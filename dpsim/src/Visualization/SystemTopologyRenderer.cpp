// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim/SystemTopologyRenderer.h>

#include <dpsim-models/Logger.h>
#include <dpsim-models/TopologicalPowerComp.h>

#ifdef WITH_GRAPHVIZ
#include <dpsim-models/Graph.h>
#endif

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <queue>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace CPS;

namespace CPS::TopologyRenderDetail {

enum class SymbolOrientation { Horizontal, Vertical };

enum class BusbarOrientation { Vertical, Horizontal };

struct SvgBounds {
  double minX = 0.0;
  double minY = 0.0;
  double maxX = 0.0;
  double maxY = 0.0;
  Bool valid = false;

  double width() const { return maxX - minX; }

  double height() const { return maxY - minY; }
};

struct Direction {
  double x = 0.0;
  double y = 0.0;
};

struct LayoutResult {
  Bool customLayoutApplied = false;
  Bool circular = false;

  std::unordered_set<String> backboneBuses;
  std::unordered_set<String> backboneComponents;

  // Direction of the main electrical path through a bus.
  std::unordered_map<String, Direction> mainDirections;

  // Preferred free direction for leaves / shunts / external branches.
  std::unordered_map<String, Direction> leafDirections;
};

struct RenderStyle {
  static constexpr Bool showNames = true;
  static constexpr Bool showTypes = false;
  static constexpr Bool useSymbols = true;
  static constexpr Bool collapseLines = false;

  static constexpr Real componentWidth = 0.60;
  static constexpr Real componentHeight = 0.60;

  static constexpr Real busWidth = 0.10;
  static constexpr Real busHeight = 0.65;

  static constexpr Real nodeSeparation = 0.44;
  static constexpr Real rankSeparation = 0.42;
};

struct GroundPlacement {
  std::filesystem::path symbol;
  SvgBounds bounds;

  double lineX0 = 0.0;
  double lineY0 = 0.0;
  double lineX1 = 0.0;
  double lineY1 = 0.0;
};

double centerX(const SvgBounds &bounds);

double centerY(const SvgBounds &bounds);

Bool extractNodeBounds(const String &svg, const String &graphNodeName,
                       SvgBounds &bounds);

Bool replaceNodeBounds(String &svg, const String &graphNodeName,
                       const SvgBounds &bounds);

Bool moveNodeCenter(String &svg, const String &graphNodeName, double x,
                    double y);

LayoutResult applyTopologyLayout(String &svg, const SystemTopology &system,
                                 SystemTopologyRenderer::Layout layout,
                                 Bool collapseLines);

Bool isLineLike(TopologicalPowerComp &component);

Bool collisionBoundsOverlap(const SvgBounds &first, const SvgBounds &second,
                            const double clearance = 0.0) {

  return first.valid && second.valid && first.maxX + clearance > second.minX &&
         first.minX - clearance < second.maxX &&
         first.maxY + clearance > second.minY &&
         first.minY - clearance < second.maxY;
}

struct LayoutCollisionObject {
  String name;
  SvgBounds bounds;
  Bool bus = false;
};

std::vector<LayoutCollisionObject>
collectLayoutCollisionObjects(const String &svg, const SystemTopology &system) {

  std::vector<LayoutCollisionObject> objects;

  for (const auto &node : system.mNodes) {

    if (!node || node->isGround()) {

      continue;
    }

    LayoutCollisionObject object;

    object.name = "bus_" + node->uid();

    object.bus = true;

    if (extractNodeBounds(svg, object.name, object.bounds)) {

      objects.push_back(object);
    }
  }

  for (const auto &identifiedObject : system.mComponents) {

    if (!identifiedObject) {
      continue;
    }

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component) {
      continue;
    }

    if (RenderStyle::collapseLines && isLineLike(*component)) {

      continue;
    }

    LayoutCollisionObject object;

    object.name = "component_" + component->uid();

    if (extractNodeBounds(svg, object.name, object.bounds)) {

      objects.push_back(object);
    }
  }

  return objects;
}

UInt countObjectOverlaps(const std::vector<LayoutCollisionObject> &objects,
                         const double clearance) {

  UInt count = 0;

  for (std::size_t first = 0; first < objects.size(); ++first) {

    for (std::size_t second = first + 1; second < objects.size(); ++second) {

      if (collisionBoundsOverlap(objects[first].bounds, objects[second].bounds,
                                 clearance)) {

        ++count;
      }
    }
  }

  return count;
}

void moveCollisionObject(LayoutCollisionObject &object, const double dx,
                         const double dy) {

  object.bounds.minX += dx;
  object.bounds.maxX += dx;
  object.bounds.minY += dy;
  object.bounds.maxY += dy;
}

void writeCollisionObjectsToSvg(
    String &svg, const std::vector<LayoutCollisionObject> &objects) {

  for (const auto &object : objects) {
    replaceNodeBounds(svg, object.name, object.bounds);
  }
}

UInt resolveObjectOverlaps(String &svg, const SystemTopology &system) {

  constexpr double clearance = 12.0;
  constexpr UInt maximumIterations = 20;

  // Parse SVG geometry once. Iterate entirely on cached bounds, then commit
  // the final positions to SVG once.
  auto objects = collectLayoutCollisionObjects(svg, system);

  if (objects.empty())
    return 0;

  UInt lastOverlapCount = std::numeric_limits<UInt>::max();

  for (UInt iteration = 0; iteration < maximumIterations; ++iteration) {

    const UInt overlapCount = countObjectOverlaps(objects, clearance);

    if (overlapCount == 0) {
      writeCollisionObjectsToSvg(svg, objects);
      return 0;
    }

    if (overlapCount >= 4 && iteration < 10) {
      double centroidX = 0.0;
      double centroidY = 0.0;

      for (const auto &object : objects) {
        centroidX += centerX(object.bounds);
        centroidY += centerY(object.bounds);
      }

      centroidX /= static_cast<double>(objects.size());
      centroidY /= static_cast<double>(objects.size());

      const double scale = overlapCount > 15 ? 1.14 : 1.09;

      for (auto &object : objects) {
        const double x = centerX(object.bounds);
        const double y = centerY(object.bounds);
        moveCollisionObject(object, (x - centroidX) * (scale - 1.0),
                            (y - centroidY) * (scale - 1.0));
      }
    }

    Bool movedSomething = false;
    Bool busBusOverlap = false;

    for (std::size_t first = 0; first < objects.size(); ++first) {
      for (std::size_t second = first + 1; second < objects.size(); ++second) {

        auto &firstObject = objects[first];
        auto &secondObject = objects[second];

        if (!collisionBoundsOverlap(firstObject.bounds, secondObject.bounds,
                                    clearance))
          continue;

        const double requiredX =
            std::min(firstObject.bounds.maxX, secondObject.bounds.maxX) -
            std::max(firstObject.bounds.minX, secondObject.bounds.minX) +
            2.0 * clearance;

        const double requiredY =
            std::min(firstObject.bounds.maxY, secondObject.bounds.maxY) -
            std::max(firstObject.bounds.minY, secondObject.bounds.minY) +
            2.0 * clearance;

        const Bool separateHorizontally = requiredX <= requiredY;
        const double firstCenter = separateHorizontally
                                       ? centerX(firstObject.bounds)
                                       : centerY(firstObject.bounds);
        const double secondCenter = separateHorizontally
                                        ? centerX(secondObject.bounds)
                                        : centerY(secondObject.bounds);

        double direction = secondCenter > firstCenter ? 1.0 : -1.0;
        if (std::abs(secondCenter - firstCenter) < 1.0e-7)
          direction = firstObject.name < secondObject.name ? 1.0 : -1.0;

        const double displacement =
            separateHorizontally ? requiredX : requiredY;

        if (firstObject.bus && secondObject.bus) {
          // Preserve the topology-aware backbone; expand globally instead.
          busBusOverlap = true;
          continue;
        }

        if (firstObject.bus && !secondObject.bus) {
          moveCollisionObject(
              secondObject,
              separateHorizontally ? direction * displacement : 0.0,
              separateHorizontally ? 0.0 : direction * displacement);
        } else if (!firstObject.bus && secondObject.bus) {
          moveCollisionObject(
              firstObject,
              separateHorizontally ? -direction * displacement : 0.0,
              separateHorizontally ? 0.0 : -direction * displacement);
        } else {
          moveCollisionObject(
              firstObject,
              separateHorizontally ? -direction * 0.5 * displacement : 0.0,
              separateHorizontally ? 0.0 : -direction * 0.5 * displacement);
          moveCollisionObject(
              secondObject,
              separateHorizontally ? direction * 0.5 * displacement : 0.0,
              separateHorizontally ? 0.0 : direction * 0.5 * displacement);
        }
        movedSomething = true;
      }
    }

    if (busBusOverlap) {
      double centroidX = 0.0;
      double centroidY = 0.0;
      for (const auto &object : objects) {
        centroidX += centerX(object.bounds);
        centroidY += centerY(object.bounds);
      }
      centroidX /= static_cast<double>(objects.size());
      centroidY /= static_cast<double>(objects.size());

      constexpr double topologyScale = 1.08;
      for (auto &object : objects) {
        const double x = centerX(object.bounds);
        const double y = centerY(object.bounds);
        moveCollisionObject(object, (x - centroidX) * (topologyScale - 1.0),
                            (y - centroidY) * (topologyScale - 1.0));
      }
      movedSomething = true;
    }

    if (!movedSomething || overlapCount == lastOverlapCount) {
      double centroidX = 0.0;
      double centroidY = 0.0;
      for (const auto &object : objects) {
        centroidX += centerX(object.bounds);
        centroidY += centerY(object.bounds);
      }
      centroidX /= static_cast<double>(objects.size());
      centroidY /= static_cast<double>(objects.size());

      for (auto &object : objects) {
        const double x = centerX(object.bounds);
        const double y = centerY(object.bounds);
        moveCollisionObject(object, (x - centroidX) * 0.04,
                            (y - centroidY) * 0.04);
      }
    }

    lastOverlapCount = overlapCount;
  }

  const UInt remaining = countObjectOverlaps(objects, clearance);
  writeCollisionObjectsToSvg(svg, objects);
  return remaining;
}

void resolveGroundOverlaps(std::vector<GroundPlacement> &grounds,
                           const std::vector<SvgBounds> &obstacles) {

  constexpr double clearance = 6.0;

  constexpr double step = 8.0;

  for (std::size_t index = 0; index < grounds.size(); ++index) {

    auto &ground = grounds[index];

    double directionX = ground.lineX1 - ground.lineX0;

    double directionY = ground.lineY1 - ground.lineY0;

    const double magnitude =
        std::sqrt(directionX * directionX + directionY * directionY);

    if (magnitude > 1.0e-9) {

      directionX /= magnitude;

      directionY /= magnitude;
    }

    for (UInt iteration = 0; iteration < 20; ++iteration) {

      Bool collision = false;

      for (const auto &obstacle : obstacles) {

        if (collisionBoundsOverlap(ground.bounds, obstacle, clearance)) {

          collision = true;

          break;
        }
      }

      if (!collision) {

        for (std::size_t previous = 0; previous < index; ++previous) {

          if (collisionBoundsOverlap(ground.bounds, grounds[previous].bounds,
                                     clearance)) {

            collision = true;

            break;
          }
        }
      }

      if (!collision) {
        break;
      }

      const double dx = directionX * step;

      const double dy = directionY * step;

      ground.bounds.minX += dx;

      ground.bounds.maxX += dx;

      ground.bounds.minY += dy;

      ground.bounds.maxY += dy;

      ground.lineX1 += dx;

      ground.lineY1 += dy;
    }
  }
}

String postProcessSvg(String svg, const SystemTopology &system,
                      const std::filesystem::path &symbolDirectory,
                      const SystemTopologyRenderer::Options &options,
                      const LayoutResult &layoutResult);

} // namespace CPS::TopologyRenderDetail

// =============================================================================
// Topology analysis and placement
// =============================================================================

namespace CPS::TopologyRenderDetail {
namespace {

struct SeriesEdge {
  const TopologicalNode *first = nullptr;
  const TopologicalNode *second = nullptr;
  String componentNodeName;
};

struct Model {
  std::vector<const TopologicalNode *> buses;

  std::unordered_map<const TopologicalNode *, String> busNodeNames;

  std::unordered_map<const TopologicalNode *, std::vector<std::size_t>>
      incidentSeriesEdges;

  std::vector<SeriesEdge> seriesEdges;

  std::unordered_map<const TopologicalNode *, std::vector<String>>
      terminalComponents;
};

Direction layoutNormalize(Direction direction) {
  const double magnitude =
      std::sqrt(direction.x * direction.x + direction.y * direction.y);

  if (magnitude <= 1.0e-12)
    return {1.0, 0.0};

  direction.x /= magnitude;
  direction.y /= magnitude;
  return direction;
}

Direction layoutCardinal(Direction direction) {
  direction = layoutNormalize(direction);

  if (std::abs(direction.x) >= std::abs(direction.y))
    return {direction.x >= 0.0 ? 1.0 : -1.0, 0.0};

  return {0.0, direction.y >= 0.0 ? 1.0 : -1.0};
}

Direction layoutPerpendicular(Direction direction) {
  direction = layoutNormalize(direction);
  return {-direction.y, direction.x};
}

const TopologicalNode *otherEnd(const SeriesEdge &edge,
                                const TopologicalNode *node) {

  if (edge.first == node)
    return edge.second;

  if (edge.second == node)
    return edge.first;

  return nullptr;
}

Model buildModel(const SystemTopology &system) {
  Model model;

  for (const auto &node : system.mNodes) {
    if (!node || node->isGround())
      continue;

    model.buses.push_back(node.get());
    model.busNodeNames[node.get()] = "bus_" + node->uid();
    model.incidentSeriesEdges[node.get()];
  }

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    std::vector<const TopologicalNode *> nonGroundNodes;
    Bool hasGround = false;

    for (const auto &terminal : component->topologicalTerminals()) {
      if (!terminal)
        continue;

      const auto node = terminal->topologicalNodes();

      if (!node)
        continue;

      if (node->isGround()) {
        hasGround = true;
        continue;
      }

      nonGroundNodes.push_back(node.get());
    }

    const String componentNodeName = "component_" + component->uid();

    if (nonGroundNodes.size() == 2 && !hasGround) {
      const std::size_t edgeIndex = model.seriesEdges.size();

      model.seriesEdges.push_back(
          {nonGroundNodes[0], nonGroundNodes[1], componentNodeName});

      model.incidentSeriesEdges[nonGroundNodes[0]].push_back(edgeIndex);

      model.incidentSeriesEdges[nonGroundNodes[1]].push_back(edgeIndex);

    } else if (nonGroundNodes.size() == 1 && !hasGround) {
      model.terminalComponents[nonGroundNodes.front()].push_back(
          componentNodeName);
    }
  }

  return model;
}

String componentBetween(const Model &model, const TopologicalNode *first,
                        const TopologicalNode *second) {

  const auto iterator = model.incidentSeriesEdges.find(first);

  if (iterator == model.incidentSeriesEdges.end())
    return {};

  for (const std::size_t edgeIndex : iterator->second) {
    const auto &edge = model.seriesEdges[edgeIndex];

    if (otherEnd(edge, first) == second)
      return edge.componentNodeName;
  }

  return {};
}

std::vector<const TopologicalNode *> neighbors(const Model &model,
                                               const TopologicalNode *node) {

  std::vector<const TopologicalNode *> result;

  const auto iterator = model.incidentSeriesEdges.find(node);

  if (iterator == model.incidentSeriesEdges.end())
    return result;

  result.reserve(iterator->second.size());

  for (const std::size_t edgeIndex : iterator->second) {
    const auto *other = otherEnd(model.seriesEdges[edgeIndex], node);

    if (other)
      result.push_back(other);
  }

  return result;
}

// -----------------------------------------------------------------------------
// General meshed-core detection
//
// Recursive degree-0/1 pruning leaves the graph's 2-core. Unlike the earlier
// simple-cycle check, the 2-core may contain many interconnected cycles. This is
// the correct electrical backbone for systems such as IEEE-39.
// -----------------------------------------------------------------------------

std::unordered_set<const TopologicalNode *> detectTwoCore(const Model &model) {

  std::unordered_map<const TopologicalNode *, UInt> degree;

  std::queue<const TopologicalNode *> queue;

  std::unordered_set<const TopologicalNode *> removed;

  for (const auto *bus : model.buses) {

    degree[bus] = static_cast<UInt>(model.incidentSeriesEdges.at(bus).size());

    if (degree[bus] < 2) {

      queue.push(bus);
    }
  }

  while (!queue.empty()) {

    const auto *bus = queue.front();

    queue.pop();

    if (removed.find(bus) != removed.end()) {

      continue;
    }

    removed.insert(bus);

    for (const auto *neighbor : neighbors(model, bus)) {

      if (removed.find(neighbor) != removed.end()) {

        continue;
      }

      auto &neighborDegree = degree[neighbor];

      if (neighborDegree > 0) {

        --neighborDegree;
      }

      if (neighborDegree == 1) {

        queue.push(neighbor);
      }
    }
  }

  std::unordered_set<const TopologicalNode *> core;

  for (const auto *bus : model.buses) {

    if (removed.find(bus) == removed.end()) {

      core.insert(bus);
    }
  }

  return core;
}

Bool isSimpleCycleCore(
    const Model &model,
    const std::unordered_set<const TopologicalNode *> &core) {

  if (core.size() < 3) {

    return false;
  }

  for (const auto *bus : core) {

    UInt coreDegree = 0;

    for (const auto *neighbor : neighbors(model, bus)) {

      if (core.find(neighbor) != core.end()) {

        ++coreDegree;
      }
    }

    if (coreDegree != 2) {

      return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Biconnected decomposition of the meshed 2-core.
//
// The renderer does not place each electrical bus independently anymore.
// Biconnected blocks are detected explicitly and used by the force layout as
// cohesive topological units. Articulation buses naturally tie neighbouring
// blocks together.
// -----------------------------------------------------------------------------

struct CoreBlock {
  std::vector<const TopologicalNode *> buses;
};

std::vector<CoreBlock>
biconnectedCoreBlocks(const Model &model,
                      const std::unordered_set<const TopologicalNode *> &core) {

  struct DfsEdge {
    const TopologicalNode *first = nullptr;

    const TopologicalNode *second = nullptr;
  };

  std::unordered_map<const TopologicalNode *, Int> discovery;

  std::unordered_map<const TopologicalNode *, Int> low;

  std::unordered_map<const TopologicalNode *, const TopologicalNode *> parent;

  std::vector<DfsEdge> edgeStack;

  std::vector<CoreBlock> blocks;

  Int time = 0;

  std::function<void(const TopologicalNode *)> dfs;

  auto emitBlock = [&](const TopologicalNode *stopFirst,
                       const TopologicalNode *stopSecond) {
    std::unordered_set<const TopologicalNode *> members;

    while (!edgeStack.empty()) {

      const DfsEdge edge = edgeStack.back();

      edgeStack.pop_back();

      members.insert(edge.first);

      members.insert(edge.second);

      if ((edge.first == stopFirst && edge.second == stopSecond) ||
          (edge.first == stopSecond && edge.second == stopFirst)) {

        break;
      }
    }

    if (members.size() >= 2) {

      CoreBlock block;

      block.buses.assign(members.begin(), members.end());

      blocks.push_back(std::move(block));
    }
  };

  dfs = [&](const TopologicalNode *bus) {
    discovery[bus] = ++time;

    low[bus] = discovery[bus];

    for (const auto *neighbor : neighbors(model, bus)) {

      if (core.find(neighbor) == core.end()) {

        continue;
      }

      if (discovery.find(neighbor) == discovery.end()) {

        parent[neighbor] = bus;

        edgeStack.push_back({bus, neighbor});

        dfs(neighbor);

        low[bus] = std::min(low[bus], low[neighbor]);

        if (low[neighbor] >= discovery[bus]) {

          emitBlock(bus, neighbor);
        }

      } else if (parent[bus] != neighbor &&
                 discovery[neighbor] < discovery[bus]) {

        low[bus] = std::min(low[bus], discovery[neighbor]);

        edgeStack.push_back({bus, neighbor});
      }
    }
  };

  for (const auto *bus : core) {

    if (discovery.find(bus) != discovery.end()) {

      continue;
    }

    parent[bus] = nullptr;

    dfs(bus);

    if (!edgeStack.empty()) {

      std::unordered_set<const TopologicalNode *> members;

      while (!edgeStack.empty()) {

        members.insert(edgeStack.back().first);

        members.insert(edgeStack.back().second);

        edgeStack.pop_back();
      }

      if (members.size() >= 2) {

        CoreBlock block;

        block.buses.assign(members.begin(), members.end());

        blocks.push_back(std::move(block));
      }
    }
  }

  return blocks;
}

std::vector<const TopologicalNode *>
coreDiameterPath(const Model &model,
                 const std::unordered_set<const TopologicalNode *> &core) {

  std::vector<const TopologicalNode *> bestPath;

  for (const auto *start : core) {

    std::queue<const TopologicalNode *> queue;

    std::unordered_map<const TopologicalNode *, Int> distance;

    std::unordered_map<const TopologicalNode *, const TopologicalNode *>
        predecessor;

    distance[start] = 0;

    predecessor[start] = nullptr;

    queue.push(start);

    while (!queue.empty()) {

      const auto *bus = queue.front();

      queue.pop();

      for (const auto *neighbor : neighbors(model, bus)) {

        if (core.find(neighbor) == core.end() ||
            distance.find(neighbor) != distance.end()) {

          continue;
        }

        distance[neighbor] = distance[bus] + 1;

        predecessor[neighbor] = bus;

        queue.push(neighbor);
      }
    }

    const TopologicalNode *farthest = start;

    Int farthestDistance = 0;

    for (const auto &entry : distance) {

      if (entry.second > farthestDistance) {

        farthestDistance = entry.second;

        farthest = entry.first;
      }
    }

    std::vector<const TopologicalNode *> path;

    for (const TopologicalNode *bus = farthest; bus != nullptr;
         bus = predecessor[bus]) {

      path.push_back(bus);
    }

    std::reverse(path.begin(), path.end());

    if (path.size() > bestPath.size()) {

      bestPath = std::move(path);
    }
  }

  return bestPath;
}

std::unordered_map<const TopologicalNode *, Int>
coreDistances(const Model &model,
              const std::unordered_set<const TopologicalNode *> &core,
              const TopologicalNode *start) {

  std::unordered_map<const TopologicalNode *, Int> distance;

  if (!start) {
    return distance;
  }

  std::queue<const TopologicalNode *> queue;

  distance[start] = 0;

  queue.push(start);

  while (!queue.empty()) {

    const auto *bus = queue.front();

    queue.pop();

    for (const auto *neighbor : neighbors(model, bus)) {

      if (core.find(neighbor) == core.end() ||
          distance.find(neighbor) != distance.end()) {

        continue;
      }

      distance[neighbor] = distance[bus] + 1;

      queue.push(neighbor);
    }
  }

  return distance;
}

// -----------------------------------------------------------------------------
// Simple cycle detection
//
// Repeatedly prune degree-0/1 buses. If the remaining connected 2-core has
// degree exactly two at every bus, it is one simple closed ring.
// -----------------------------------------------------------------------------

std::vector<const TopologicalNode *> detectSimpleCycle(const Model &model) {

  if (model.buses.size() < 3)
    return {};

  std::unordered_map<const TopologicalNode *, UInt> degree;
  std::queue<const TopologicalNode *> queue;
  std::unordered_set<const TopologicalNode *> removed;

  for (const auto *bus : model.buses) {
    degree[bus] = static_cast<UInt>(model.incidentSeriesEdges.at(bus).size());

    if (degree[bus] < 2)
      queue.push(bus);
  }

  while (!queue.empty()) {
    const auto *bus = queue.front();
    queue.pop();

    if (removed.find(bus) != removed.end())
      continue;

    removed.insert(bus);

    for (const auto *neighbor : neighbors(model, bus)) {
      if (removed.find(neighbor) != removed.end())
        continue;

      auto &neighborDegree = degree[neighbor];

      if (neighborDegree > 0)
        --neighborDegree;

      if (neighborDegree == 1)
        queue.push(neighbor);
    }
  }

  std::unordered_set<const TopologicalNode *> core;

  for (const auto *bus : model.buses) {
    if (removed.find(bus) == removed.end())
      core.insert(bus);
  }

  if (core.size() < 3)
    return {};

  for (const auto *bus : core) {
    UInt coreDegree = 0;

    for (const auto *neighbor : neighbors(model, bus)) {
      if (core.find(neighbor) != core.end())
        ++coreDegree;
    }

    if (coreDegree != 2)
      return {};
  }

  const TopologicalNode *start = nullptr;

  for (const auto *bus : model.buses) {
    if (core.find(bus) != core.end()) {
      start = bus;
      break;
    }
  }

  if (!start)
    return {};

  std::vector<const TopologicalNode *> ordered;
  const TopologicalNode *previous = nullptr;
  const TopologicalNode *current = start;

  while (true) {
    ordered.push_back(current);

    const TopologicalNode *next = nullptr;

    for (const auto *neighbor : neighbors(model, current)) {
      if (core.find(neighbor) == core.end() || neighbor == previous)
        continue;

      next = neighbor;
      break;
    }

    if (!next)
      return {};

    if (next == start)
      break;

    previous = current;
    current = next;

    if (ordered.size() > core.size())
      return {};
  }

  if (ordered.size() != core.size())
    return {};

  return ordered;
}

// -----------------------------------------------------------------------------
// Longest connection / diameter
//
// For a tree this is the exact longest bus-to-bus path.
// For a non-circular meshed graph it is the longest shortest path, which is a
// deterministic and useful publication backbone without solving the NP-hard
// longest-simple-path problem.
// -----------------------------------------------------------------------------

std::vector<const TopologicalNode *> longestConnection(const Model &model) {

  if (model.buses.empty())
    return {};

  std::vector<const TopologicalNode *> bestPath;

  for (const auto *start : model.buses) {
    std::queue<const TopologicalNode *> queue;

    std::unordered_map<const TopologicalNode *, Int> distance;
    std::unordered_map<const TopologicalNode *, const TopologicalNode *>
        predecessor;

    distance[start] = 0;
    predecessor[start] = nullptr;
    queue.push(start);

    while (!queue.empty()) {
      const auto *bus = queue.front();
      queue.pop();

      for (const auto *neighbor : neighbors(model, bus)) {
        if (distance.find(neighbor) != distance.end())
          continue;

        distance[neighbor] = distance[bus] + 1;

        predecessor[neighbor] = bus;
        queue.push(neighbor);
      }
    }

    const TopologicalNode *farthest = start;
    Int farthestDistance = 0;

    for (const auto &entry : distance) {
      if (entry.second > farthestDistance) {
        farthestDistance = entry.second;
        farthest = entry.first;
      }
    }

    std::vector<const TopologicalNode *> path;

    for (const TopologicalNode *bus = farthest; bus != nullptr;
         bus = predecessor[bus]) {

      path.push_back(bus);
    }

    std::reverse(path.begin(), path.end());

    if (path.size() > bestPath.size())
      bestPath = std::move(path);
  }

  return bestPath;
}

std::pair<double, double> graphCenter(const String &svg, const Model &model) {

  double x = 0.0;
  double y = 0.0;
  std::size_t count = 0;

  for (const auto *bus : model.buses) {
    SvgBounds bounds;

    if (!extractNodeBounds(svg, model.busNodeNames.at(bus), bounds))
      continue;

    x += centerX(bounds);
    y += centerY(bounds);
    ++count;
  }

  if (count == 0)
    return {0.0, 0.0};

  return {x / static_cast<double>(count), y / static_cast<double>(count)};
}

void placeSeriesComponentsAtMidpoints(String &svg, const Model &model) {

  for (const auto &edge : model.seriesEdges) {
    SvgBounds firstBounds;
    SvgBounds secondBounds;

    if (!extractNodeBounds(svg, model.busNodeNames.at(edge.first),
                           firstBounds) ||
        !extractNodeBounds(svg, model.busNodeNames.at(edge.second),
                           secondBounds))
      continue;

    moveNodeCenter(svg, edge.componentNodeName,
                   0.5 * (centerX(firstBounds) + centerX(secondBounds)),
                   0.5 * (centerY(firstBounds) + centerY(secondBounds)));
  }
}

void placeTerminalComponents(String &svg, const Model &model,
                             const TopologicalNode *bus, Direction direction,
                             double distance, double siblingSpacing = 34.0) {

  const auto iterator = model.terminalComponents.find(bus);

  if (iterator == model.terminalComponents.end())
    return;

  SvgBounds busBounds;

  if (!extractNodeBounds(svg, model.busNodeNames.at(bus), busBounds))
    return;

  direction = layoutCardinal(direction);
  const Direction tangent = layoutPerpendicular(direction);

  for (std::size_t index = 0; index < iterator->second.size(); ++index) {

    const double centeredIndex =
        static_cast<double>(index) -
        0.5 * static_cast<double>(iterator->second.size() - 1);

    moveNodeCenter(svg, iterator->second[index],
                   centerX(busBounds) + direction.x * distance +
                       tangent.x * centeredIndex * siblingSpacing,
                   centerY(busBounds) + direction.y * distance +
                       tangent.y * centeredIndex * siblingSpacing);
  }
}

void placeExternalTree(
    String &svg, const Model &model, const TopologicalNode *bus,
    const TopologicalNode *parent,
    const std::unordered_set<const TopologicalNode *> &blocked,
    Direction outward, LayoutResult &result,
    std::unordered_set<const TopologicalNode *> &visited, UInt depth) {

  outward = layoutCardinal(outward);

  const String busName = model.busNodeNames.at(bus);

  result.leafDirections[busName] = outward;

  result.mainDirections[busName] = outward;

  SvgBounds busBounds;

  if (!extractNodeBounds(svg, busName, busBounds))
    return;

  std::vector<const TopologicalNode *> children;

  for (const auto *neighbor : neighbors(model, bus)) {
    if (neighbor == parent || blocked.find(neighbor) != blocked.end() ||
        visited.find(neighbor) != visited.end())
      continue;

    children.push_back(neighbor);
  }

  constexpr double componentStep = 50.0;
  constexpr double busStep = 100.0;
  constexpr double siblingStep = 38.0;

  const Direction side = layoutPerpendicular(outward);

  for (std::size_t index = 0; index < children.size(); ++index) {

    const double centeredIndex = static_cast<double>(index) -
                                 0.5 * static_cast<double>(children.size() - 1);

    const double sideOffset = centeredIndex * siblingStep;

    const auto *child = children[index];

    moveNodeCenter(
        svg, componentBetween(model, bus, child),
        centerX(busBounds) + outward.x * componentStep + side.x * sideOffset,
        centerY(busBounds) + outward.y * componentStep + side.y * sideOffset);

    moveNodeCenter(
        svg, model.busNodeNames.at(child),
        centerX(busBounds) + outward.x * busStep + side.x * sideOffset,
        centerY(busBounds) + outward.y * busStep + side.y * sideOffset);

    visited.insert(child);

    placeExternalTree(svg, model, child, bus, blocked, outward, result, visited,
                      depth + 1);
  }

  if (children.empty()) {
    // Terminal device continues the branch naturally outward.
    placeTerminalComponents(svg, model, bus, outward, 54.0);

    result.leafDirections[busName] = outward;

  } else {
    // Shunts at a branching bus should leave perpendicular to the local branch.
    result.leafDirections[busName] =
        layoutCardinal(layoutPerpendicular(outward));
  }
}

LayoutResult placeCircular(String &svg, const Model &model,
                           const std::vector<const TopologicalNode *> &cycle,
                           SystemTopologyRenderer::Layout layout) {

  LayoutResult result;
  result.customLayoutApplied = true;
  result.circular = true;

  const auto [centerXValue, centerYValue] = graphCenter(svg, model);

  double radiusX = 120.0;
  double radiusY = 120.0;

  if (layout == SystemTopologyRenderer::Layout::LeftToRight) {
    radiusX = 155.0;
    radiusY = 105.0;

  } else if (layout == SystemTopologyRenderer::Layout::TopToBottom) {

    radiusX = 110.0;
    radiusY = 150.0;
  }

  constexpr double pi = 3.14159265358979323846;

  std::unordered_map<const TopologicalNode *, std::pair<double, double>>
      positions;

  for (std::size_t index = 0; index < cycle.size(); ++index) {

    const double angle = -0.5 * pi + 2.0 * pi * static_cast<double>(index) /
                                         static_cast<double>(cycle.size());

    const double x = centerXValue + radiusX * std::cos(angle);

    const double y = centerYValue + radiusY * std::sin(angle);

    positions[cycle[index]] = {x, y};

    moveNodeCenter(svg, model.busNodeNames.at(cycle[index]), x, y);

    result.backboneBuses.insert(model.busNodeNames.at(cycle[index]));
  }

  for (std::size_t index = 0; index < cycle.size(); ++index) {

    const auto *first = cycle[index];
    const auto *second = cycle[(index + 1) % cycle.size()];

    const String component = componentBetween(model, first, second);

    result.backboneComponents.insert(component);

    moveNodeCenter(svg, component,
                   0.5 * (positions[first].first + positions[second].first),
                   0.5 * (positions[first].second + positions[second].second));
  }

  std::unordered_set<const TopologicalNode *> core(cycle.begin(), cycle.end());

  // Main tangent and outside normal at every cycle bus.
  for (std::size_t index = 0; index < cycle.size(); ++index) {

    const auto *previous = cycle[(index + cycle.size() - 1) % cycle.size()];

    const auto *bus = cycle[index];

    const auto *next = cycle[(index + 1) % cycle.size()];

    Direction tangent{positions[next].first - positions[previous].first,
                      positions[next].second - positions[previous].second};

    tangent = layoutNormalize(tangent);

    Direction normal = layoutPerpendicular(tangent);

    const Direction radial{positions[bus].first - centerXValue,
                           positions[bus].second - centerYValue};

    if (normal.x * radial.x + normal.y * radial.y < 0.0) {

      normal.x *= -1.0;
      normal.y *= -1.0;
    }

    const String busName = model.busNodeNames.at(bus);

    result.mainDirections[busName] = tangent;

    result.leafDirections[busName] = layoutCardinal(normal);
  }

  // External trees and terminal devices are placed outside the ring.
  std::unordered_set<const TopologicalNode *> visited;

  for (const auto *cycleBus : cycle) {
    const String busName = model.busNodeNames.at(cycleBus);

    const Direction outward = result.leafDirections.at(busName);

    SvgBounds busBounds;

    if (!extractNodeBounds(svg, busName, busBounds))
      continue;

    std::vector<const TopologicalNode *> children;

    for (const auto *neighbor : neighbors(model, cycleBus)) {
      if (core.find(neighbor) == core.end())
        children.push_back(neighbor);
    }

    const Direction side = layoutPerpendicular(outward);

    constexpr double componentStep = 50.0;
    constexpr double busStep = 100.0;
    constexpr double siblingStep = 38.0;

    for (std::size_t index = 0; index < children.size(); ++index) {

      const double centeredIndex =
          static_cast<double>(index) -
          0.5 * static_cast<double>(children.size() - 1);

      const double offset = centeredIndex * siblingStep;

      const auto *child = children[index];

      moveNodeCenter(
          svg, componentBetween(model, cycleBus, child),
          centerX(busBounds) + outward.x * componentStep + side.x * offset,
          centerY(busBounds) + outward.y * componentStep + side.y * offset);

      moveNodeCenter(svg, model.busNodeNames.at(child),
                     centerX(busBounds) + outward.x * busStep + side.x * offset,
                     centerY(busBounds) + outward.y * busStep +
                         side.y * offset);

      visited.insert(child);

      placeExternalTree(svg, model, child, cycleBus, core, outward, result,
                        visited, 1);
    }

    if (children.empty()) {
      placeTerminalComponents(svg, model, cycleBus, outward, 54.0);
    }
  }

  return result;
}

LayoutResult
placeMeshedCore(String &svg, const Model &model,
                const std::unordered_set<const TopologicalNode *> &core,
                const SystemTopologyRenderer::Layout layout) {

  LayoutResult result;

  if (core.empty()) {
    return result;
  }

  result.customLayoutApplied = true;

  result.circular = false;

  const auto blocks = biconnectedCoreBlocks(model, core);

  const auto diameter = coreDiameterPath(model, core);

  const TopologicalNode *diameterStart =
      diameter.empty() ? *core.begin() : diameter.front();

  const auto distanceFromStart = coreDistances(model, core, diameterStart);

  std::unordered_map<const TopologicalNode *, std::pair<double, double>>
      position;

  std::unordered_map<const TopologicalNode *, std::pair<double, double>>
      initialPosition;

  double initialCenterX = 0.0;

  double initialCenterY = 0.0;

  UInt initialCount = 0;

  for (const auto *bus : core) {

    SvgBounds bounds;

    if (!extractNodeBounds(svg, model.busNodeNames.at(bus), bounds)) {

      continue;
    }

    const double x = centerX(bounds);

    const double y = centerY(bounds);

    initialPosition[bus] = {x, y};

    initialCenterX += x;

    initialCenterY += y;

    ++initialCount;
  }

  if (initialCount > 0) {

    initialCenterX /= static_cast<double>(initialCount);

    initialCenterY /= static_cast<double>(initialCount);
  }

  const Bool leftToRight =
      layout == SystemTopologyRenderer::Layout::LeftToRight;

  const Bool topToBottom =
      layout == SystemTopologyRenderer::Layout::TopToBottom;

  constexpr double primarySpacing = 118.0;

  constexpr double initialSecondaryScale = 0.65;

  // ---------------------------------------------------------------------------
  // Initial block-aware embedding.
  //
  // LR/TB use the diameter only to establish the principal AXIS. The complete
  // meshed core is still laid out by a force system; it is no longer forced
  // onto one longest path.
  //
  // Quadratic starts from Graphviz geometry and is pulled toward a square.
  // ---------------------------------------------------------------------------

  for (const auto *bus : core) {

    const auto initialIterator = initialPosition.find(bus);

    const double relativeX =
        initialIterator == initialPosition.end()
            ? 0.0
            : initialIterator->second.first - initialCenterX;

    const double relativeY =
        initialIterator == initialPosition.end()
            ? 0.0
            : initialIterator->second.second - initialCenterY;

    if (leftToRight || topToBottom) {

      const auto distanceIterator = distanceFromStart.find(bus);

      const double level = distanceIterator == distanceFromStart.end()
                               ? 0.0
                               : static_cast<double>(distanceIterator->second);

      if (leftToRight) {

        position[bus] = {level * primarySpacing,
                         relativeY * initialSecondaryScale};

      } else {

        position[bus] = {relativeX * initialSecondaryScale,
                         level * primarySpacing};
      }

    } else {

      position[bus] = {relativeX * 0.85, relativeY * 0.85};
    }
  }

  struct CoreEdge {
    const TopologicalNode *first = nullptr;

    const TopologicalNode *second = nullptr;
  };

  std::vector<CoreEdge> coreEdges;

  for (const auto &edge : model.seriesEdges) {

    if (core.find(edge.first) != core.end() &&
        core.find(edge.second) != core.end()) {

      coreEdges.push_back({edge.first, edge.second});
    }
  }

  // Number of biconnected blocks touching each bus. Articulation buses are
  // deliberately kept more stable by the force solver.
  std::unordered_map<const TopologicalNode *, UInt> blockMultiplicity;

  for (const auto &block : blocks) {

    for (const auto *bus : block.buses) {

      ++blockMultiplicity[bus];
    }
  }

  constexpr double repulsionStrength = 15500.0;

  constexpr double springRestLength = 108.0;

  constexpr double springStrength = 0.045;

  constexpr double blockCohesion = 0.010;

  constexpr double axisStrength = 0.11;

  constexpr double gravityStrength = 0.010;

  constexpr double maximumStep = 11.0;

  constexpr UInt maximumForceIterations = 100;

  constexpr double forceConvergenceTolerance = 0.10;

  constexpr UInt forceStableIterationsRequired = 4;

  UInt forceStableIterations = 0;

  for (UInt iteration = 0; iteration < maximumForceIterations; ++iteration) {

    std::unordered_map<const TopologicalNode *, Direction> force;

    for (const auto *bus : core) {

      force[bus] = {0.0, 0.0};
    }

    // Pairwise core-bus repulsion.
    std::vector<const TopologicalNode *> coreBuses(core.begin(), core.end());

    for (std::size_t first = 0; first < coreBuses.size(); ++first) {

      for (std::size_t second = first + 1; second < coreBuses.size();
           ++second) {

        const auto *firstBus = coreBuses[first];

        const auto *secondBus = coreBuses[second];

        double dx = position[firstBus].first - position[secondBus].first;

        double dy = position[firstBus].second - position[secondBus].second;

        double distanceSquared = dx * dx + dy * dy;

        if (distanceSquared < 25.0) {

          dx += first < second ? -5.0 : 5.0;

          dy += first < second ? 5.0 : -5.0;

          distanceSquared = dx * dx + dy * dy;
        }

        const double distance = std::sqrt(distanceSquared);

        const double magnitude = repulsionStrength / distanceSquared;

        const double fx = magnitude * dx / distance;

        const double fy = magnitude * dy / distance;

        force[firstBus].x += fx;

        force[firstBus].y += fy;

        force[secondBus].x -= fx;

        force[secondBus].y -= fy;
      }
    }

    // Electrical edges are springs.
    for (const auto &edge : coreEdges) {

      double dx = position[edge.second].first - position[edge.first].first;

      double dy = position[edge.second].second - position[edge.first].second;

      const double distance = std::max(1.0, std::sqrt(dx * dx + dy * dy));

      const double magnitude = springStrength * (distance - springRestLength);

      dx /= distance;

      dy /= distance;

      force[edge.first].x += magnitude * dx;

      force[edge.first].y += magnitude * dy;

      force[edge.second].x -= magnitude * dx;

      force[edge.second].y -= magnitude * dy;
    }

    // Biconnected blocks act as cohesive topological groups. This is gentle:
    // edges still define the exact network geometry, while the blocks resist
    // the "exploded spaghetti" effect of independent node relaxation.
    for (const auto &block : blocks) {

      if (block.buses.empty()) {
        continue;
      }

      double blockCenterX = 0.0;

      double blockCenterY = 0.0;

      for (const auto *bus : block.buses) {

        blockCenterX += position[bus].first;

        blockCenterY += position[bus].second;
      }

      blockCenterX /= static_cast<double>(block.buses.size());

      blockCenterY /= static_cast<double>(block.buses.size());

      for (const auto *bus : block.buses) {

        const double articulationFactor =
            blockMultiplicity[bus] > 1 ? 0.45 : 1.0;

        force[bus].x += articulationFactor * blockCohesion *
                        (blockCenterX - position[bus].first);

        force[bus].y += articulationFactor * blockCohesion *
                        (blockCenterY - position[bus].second);
      }
    }

    // Principal-axis bias for LR/TB. It is a bias, not a hard rank, so loops
    // retain enough freedom to spread naturally.
    if (leftToRight || topToBottom) {

      for (const auto *bus : core) {

        const auto distanceIterator = distanceFromStart.find(bus);

        if (distanceIterator == distanceFromStart.end()) {

          continue;
        }

        const double target =
            static_cast<double>(distanceIterator->second) * primarySpacing;

        if (leftToRight) {

          force[bus].x += axisStrength * (target - position[bus].first);

        } else {

          force[bus].y += axisStrength * (target - position[bus].second);
        }
      }
    }

    // Mild gravity keeps the general meshed core compact.
    double centerXValue = 0.0;

    double centerYValue = 0.0;

    for (const auto *bus : core) {

      centerXValue += position[bus].first;

      centerYValue += position[bus].second;
    }

    centerXValue /= static_cast<double>(core.size());

    centerYValue /= static_cast<double>(core.size());

    double maximumMovement = 0.0;

    for (const auto *bus : core) {

      force[bus].x += gravityStrength * (centerXValue - position[bus].first);

      force[bus].y += gravityStrength * (centerYValue - position[bus].second);

      const double articulationDamping =
          blockMultiplicity[bus] > 1 ? 0.55 : 1.0;

      const double stepX = std::clamp(articulationDamping * force[bus].x,
                                      -maximumStep, maximumStep);

      const double stepY = std::clamp(articulationDamping * force[bus].y,
                                      -maximumStep, maximumStep);

      maximumMovement =
          std::max(maximumMovement, std::sqrt(stepX * stepX + stepY * stepY));

      position[bus].first += stepX;

      position[bus].second += stepY;
    }

    if (maximumMovement <= forceConvergenceTolerance) {

      ++forceStableIterations;

      if (forceStableIterations >= forceStableIterationsRequired) {

        break;
      }

    } else {

      forceStableIterations = 0;
    }
  }

  // ---------------------------------------------------------------------------
  // Final core-only overlap relaxation.
  // ---------------------------------------------------------------------------

  constexpr double minimumBusDistance = 82.0;

  for (UInt iteration = 0; iteration < 18; ++iteration) {

    Bool moved = false;

    std::vector<const TopologicalNode *> coreBuses(core.begin(), core.end());

    for (std::size_t first = 0; first < coreBuses.size(); ++first) {

      for (std::size_t second = first + 1; second < coreBuses.size();
           ++second) {

        const auto *firstBus = coreBuses[first];

        const auto *secondBus = coreBuses[second];

        double dx = position[secondBus].first - position[firstBus].first;

        double dy = position[secondBus].second - position[firstBus].second;

        double distance = std::sqrt(dx * dx + dy * dy);

        if (distance >= minimumBusDistance) {

          continue;
        }

        if (distance < 1.0) {

          dx = 1.0;

          dy = 1.0;

          distance = std::sqrt(2.0);
        }

        const double push = 0.5 * (minimumBusDistance - distance);

        dx /= distance;

        dy /= distance;

        position[firstBus].first -= dx * push;

        position[firstBus].second -= dy * push;

        position[secondBus].first += dx * push;

        position[secondBus].second += dy * push;

        moved = true;
      }
    }

    if (!moved) {
      break;
    }
  }

  // Re-centre the custom core near Graphviz's original graph centre.
  const auto [graphCenterX, graphCenterY] = graphCenter(svg, model);

  double coreCenterX = 0.0;

  double coreCenterY = 0.0;

  for (const auto *bus : core) {

    coreCenterX += position[bus].first;

    coreCenterY += position[bus].second;
  }

  coreCenterX /= static_cast<double>(core.size());

  coreCenterY /= static_cast<double>(core.size());

  for (const auto *bus : core) {

    position[bus].first += graphCenterX - coreCenterX;

    position[bus].second += graphCenterY - coreCenterY;

    moveNodeCenter(svg, model.busNodeNames.at(bus), position[bus].first,
                   position[bus].second);

    result.backboneBuses.insert(model.busNodeNames.at(bus));
  }

  // Place every component belonging to the meshed core directly between its
  // buses. External tree components are handled separately below.
  for (const auto &edge : model.seriesEdges) {

    if (core.find(edge.first) == core.end() ||
        core.find(edge.second) == core.end()) {

      continue;
    }

    moveNodeCenter(
        svg, edge.componentNodeName,
        0.5 * (position[edge.first].first + position[edge.second].first),
        0.5 * (position[edge.first].second + position[edge.second].second));

    result.backboneComponents.insert(edge.componentNodeName);
  }

  double placedCoreCenterX = 0.0;

  double placedCoreCenterY = 0.0;

  for (const auto *bus : core) {

    placedCoreCenterX += position[bus].first;

    placedCoreCenterY += position[bus].second;
  }

  placedCoreCenterX /= static_cast<double>(core.size());

  placedCoreCenterY /= static_cast<double>(core.size());

  // Determine a local principal direction and a preferred OUTWARD leaf sector
  // for every meshed-core bus.
  for (const auto *bus : core) {

    std::vector<const TopologicalNode *> coreNeighbors;

    for (const auto *neighbor : neighbors(model, bus)) {

      if (core.find(neighbor) != core.end()) {

        coreNeighbors.push_back(neighbor);
      }
    }

    Direction mainDirection{1.0, 0.0};

    if (coreNeighbors.size() >= 2) {

      double bestDistance = -1.0;

      for (std::size_t first = 0; first < coreNeighbors.size(); ++first) {

        for (std::size_t second = first + 1; second < coreNeighbors.size();
             ++second) {

          const double dx = position[coreNeighbors[second]].first -
                            position[coreNeighbors[first]].first;

          const double dy = position[coreNeighbors[second]].second -
                            position[coreNeighbors[first]].second;

          const double distance = dx * dx + dy * dy;

          if (distance > bestDistance) {

            bestDistance = distance;

            mainDirection = {dx, dy};
          }
        }
      }

    } else if (coreNeighbors.size() == 1) {

      mainDirection = {
          position[coreNeighbors.front()].first - position[bus].first,
          position[coreNeighbors.front()].second - position[bus].second};
    }

    mainDirection = layoutNormalize(mainDirection);

    Direction outward{position[bus].first - placedCoreCenterX,
                      position[bus].second - placedCoreCenterY};

    if (std::sqrt(outward.x * outward.x + outward.y * outward.y) < 20.0) {

      outward = layoutPerpendicular(mainDirection);
    }

    outward = layoutCardinal(outward);

    const String busName = model.busNodeNames.at(bus);

    result.mainDirections[busName] = mainDirection;

    result.leafDirections[busName] = outward;
  }

  // ---------------------------------------------------------------------------
  // Everything outside the 2-core is a radial tree. Place complete subtrees in
  // the OUTWARD sector of their attachment bus.
  // ---------------------------------------------------------------------------

  std::unordered_set<const TopologicalNode *> visited;

  for (const auto *coreBus : core) {

    const String busName = model.busNodeNames.at(coreBus);

    const Direction outward = result.leafDirections.at(busName);

    SvgBounds busBounds;

    if (!extractNodeBounds(svg, busName, busBounds)) {

      continue;
    }

    std::vector<const TopologicalNode *> externalChildren;

    for (const auto *neighbor : neighbors(model, coreBus)) {

      if (core.find(neighbor) == core.end()) {

        externalChildren.push_back(neighbor);
      }
    }

    const Direction side = layoutPerpendicular(outward);

    constexpr double componentStep = 58.0;

    constexpr double busStep = 118.0;

    constexpr double siblingStep = 46.0;

    for (std::size_t index = 0; index < externalChildren.size(); ++index) {

      const double centeredIndex =
          static_cast<double>(index) -
          0.5 * static_cast<double>(externalChildren.size() - 1);

      const double sideOffset = centeredIndex * siblingStep;

      const auto *child = externalChildren[index];

      moveNodeCenter(
          svg, componentBetween(model, coreBus, child),
          centerX(busBounds) + outward.x * componentStep + side.x * sideOffset,
          centerY(busBounds) + outward.y * componentStep + side.y * sideOffset);

      moveNodeCenter(
          svg, model.busNodeNames.at(child),
          centerX(busBounds) + outward.x * busStep + side.x * sideOffset,
          centerY(busBounds) + outward.y * busStep + side.y * sideOffset);

      visited.insert(child);

      placeExternalTree(svg, model, child, coreBus, core, outward, result,
                        visited, 1);
    }

    // Generators / injections directly attached to a core bus also occupy the
    // outward sector. Local collision relaxation can subsequently split several
    // terminal devices into neighbouring lanes.
    placeTerminalComponents(svg, model, coreBus, outward, 66.0, 44.0);
  }

  return result;
}

LayoutResult
placeLongestConnection(String &svg, const Model &model,
                       const std::vector<const TopologicalNode *> &backbone,
                       SystemTopologyRenderer::Layout layout) {

  LayoutResult result;

  if (backbone.empty())
    return result;

  result.customLayoutApplied = true;

  const auto [centerXValue, centerYValue] = graphCenter(svg, model);

  const Bool leftToRight =
      layout == SystemTopologyRenderer::Layout::LeftToRight;

  const Direction primary =
      leftToRight ? Direction{1.0, 0.0} : Direction{0.0, 1.0};

  const Direction side = layoutPerpendicular(primary);

  constexpr double busSpacing = 112.0;

  const double totalLength =
      static_cast<double>(backbone.size() - 1) * busSpacing;

  const double startX = centerXValue - 0.5 * totalLength * primary.x;

  const double startY = centerYValue - 0.5 * totalLength * primary.y;

  std::unordered_set<const TopologicalNode *> backboneSet(backbone.begin(),
                                                          backbone.end());

  for (std::size_t index = 0; index < backbone.size(); ++index) {

    const double x =
        startX + static_cast<double>(index) * busSpacing * primary.x;

    const double y =
        startY + static_cast<double>(index) * busSpacing * primary.y;

    moveNodeCenter(svg, model.busNodeNames.at(backbone[index]), x, y);

    const String busName = model.busNodeNames.at(backbone[index]);

    result.backboneBuses.insert(busName);
    result.mainDirections[busName] = primary;

    // Internal backbone buses use a perpendicular free side for shunts.
    // Endpoint buses continue naturally beyond the end.
    if (index == 0) {
      result.leafDirections[busName] = {-primary.x, -primary.y};

    } else if (index + 1 == backbone.size()) {
      result.leafDirections[busName] = primary;

    } else {
      result.leafDirections[busName] = side;
    }
  }

  // Put backbone components exactly between adjacent buses.
  for (std::size_t index = 0; index + 1 < backbone.size(); ++index) {

    const String component =
        componentBetween(model, backbone[index], backbone[index + 1]);

    result.backboneComponents.insert(component);

    SvgBounds firstBounds;
    SvgBounds secondBounds;

    if (!extractNodeBounds(svg, model.busNodeNames.at(backbone[index]),
                           firstBounds) ||
        !extractNodeBounds(svg, model.busNodeNames.at(backbone[index + 1]),
                           secondBounds))
      continue;

    moveNodeCenter(svg, component,
                   0.5 * (centerX(firstBounds) + centerX(secondBounds)),
                   0.5 * (centerY(firstBounds) + centerY(secondBounds)));
  }

  // Off-backbone trees become leaves from the relevant middle bus.
  std::unordered_set<const TopologicalNode *> visited;

  for (std::size_t backboneIndex = 0; backboneIndex < backbone.size();
       ++backboneIndex) {

    const auto *bus = backbone[backboneIndex];

    SvgBounds busBounds;

    if (!extractNodeBounds(svg, model.busNodeNames.at(bus), busBounds))
      continue;

    std::vector<const TopologicalNode *> branches;

    for (const auto *neighbor : neighbors(model, bus)) {
      if (backboneSet.find(neighbor) == backboneSet.end())
        branches.push_back(neighbor);
    }

    // Alternate branches above/below for LR, left/right for TB.
    for (std::size_t branchIndex = 0; branchIndex < branches.size();
         ++branchIndex) {

      const double sign = branchIndex % 2 == 0 ? 1.0 : -1.0;

      Direction outward{side.x * sign, side.y * sign};

      const double laneOffset = static_cast<double>(branchIndex / 2) * 34.0;

      const Direction along = primary;

      const auto *child = branches[branchIndex];

      moveNodeCenter(
          svg, componentBetween(model, bus, child),
          centerX(busBounds) + outward.x * 50.0 + along.x * laneOffset,
          centerY(busBounds) + outward.y * 50.0 + along.y * laneOffset);

      moveNodeCenter(
          svg, model.busNodeNames.at(child),
          centerX(busBounds) + outward.x * 100.0 + along.x * laneOffset,
          centerY(busBounds) + outward.y * 100.0 + along.y * laneOffset);

      visited.insert(child);

      placeExternalTree(svg, model, child, bus, backboneSet, outward, result,
                        visited, 1);
    }

    // Terminal devices on backbone endpoints continue the longest connection.
    // At middle buses they are ordinary leaves.
    if (backboneIndex == 0) {
      placeTerminalComponents(svg, model, bus, {-primary.x, -primary.y}, 54.0);

    } else if (backboneIndex + 1 == backbone.size()) {

      placeTerminalComponents(svg, model, bus, primary, 54.0);

    } else if (branches.empty()) {

      placeTerminalComponents(svg, model, bus, side, 54.0);
    }
  }

  // If a disconnected component was not moved by the selected diameter,
  // Graphviz retains its initial placement.
  placeSeriesComponentsAtMidpoints(svg, model);

  return result;
}

} // namespace

LayoutResult applyTopologyLayout(String &svg, const SystemTopology &system,
                                 const SystemTopologyRenderer::Layout layout,
                                 const Bool collapseLines) {

  // Manual positioning assumes explicit component nodes.
  if (collapseLines) {
    return {};
  }

  const Model model = buildModel(system);

  const auto core = detectTwoCore(model);

  if (!core.empty() && isSimpleCycleCore(model, core)) {

    const auto cycle = detectSimpleCycle(model);

    if (!cycle.empty()) {

      return placeCircular(svg, model, cycle, layout);
    }
  }

  if (!core.empty()) {

    return placeMeshedCore(svg, model, core, layout);
  }

  // Acyclic/radial case.
  //
  // Preserve the established compact Graphviz behaviour for small Quadratic
  // examples such as the simplified RLC renderer test. LeftToRight and
  // TopToBottom continue to use the explicit longest-connection/tree layout.
  // Larger radial Quadratic networks use the topology-aware tree algorithm.
  if (layout == SystemTopologyRenderer::Layout::Quadratic &&
      model.buses.size() <= 12) {

    return {};
  }

  return placeLongestConnection(svg, model, longestConnection(model), layout);
}

} // namespace CPS::TopologyRenderDetail

// =============================================================================
// SVG rendering, broad-side bus routing, symbols, grounds and labels
// =============================================================================

namespace CPS::TopologyRenderDetail {
namespace {

constexpr const char *kFont = "Latin Modern Roman";

constexpr double kComponentLabelSize = 5.5;
constexpr double kTypeLabelSize = 4.25;
constexpr double kBusLabelSize = 5.0;

Bool contains(const String &text, const String &needle) {

  return text.find(needle) != String::npos;
}

String lowerCopy(String value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char character) {
                   return static_cast<char>(std::tolower(character));
                 });

  return value;
}

String formatReal(const double value) {
  std::ostringstream stream;

  stream << std::fixed << std::setprecision(6) << value;

  return stream.str();
}

String readTextFile(const std::filesystem::path &path) {

  std::ifstream stream(path);

  if (!stream.is_open()) {
    throw std::runtime_error("Cannot open SVG symbol file: " + path.string());
  }

  std::ostringstream buffer;
  buffer << stream.rdbuf();

  return buffer.str();
}

void replaceAll(String &text, const String &from, const String &to) {

  if (from.empty())
    return;

  std::size_t position = 0;

  while ((position = text.find(from, position)) != String::npos) {

    text.replace(position, from.size(), to);

    position += to.size();
  }
}

void forceMonochromeSvg(String &svg) {
  replaceAll(svg, "color-scheme: light dark;", "color-scheme: light;");

  replaceAll(svg, "color-scheme:light dark;", "color-scheme:light;");

  replaceAll(svg, "light-dark(rgb(0, 0, 0), rgb(255, 255, 255))", "#000000");

  replaceAll(svg, "light-dark(rgb(0,0,0), rgb(255,255,255))", "#000000");

  replaceAll(svg, "light-dark(rgb(0, 0, 0), rgb(237, 237, 237))", "#000000");

  replaceAll(svg, "light-dark(rgb(0,0,0), rgb(237,237,237))", "#000000");

  replaceAll(svg, "light-dark(#000000, #ffffff)", "#000000");

  replaceAll(svg, "var(--ge-adaptive-bg, #ffffff)", "#ffffff");

  svg = std::regex_replace(
      svg,
      std::regex(R"regex(color-scheme\s*:\s*light\s+dark\s*;?)regex",
                 std::regex::icase),
      "color-scheme: light;");
}

SvgBounds extractPolygonBounds(const String &group) {

  SvgBounds bounds;

  std::smatch polygonMatch;

  if (!std::regex_search(group, polygonMatch,
                         std::regex(R"regex(points="([^"]+)")regex")))
    return bounds;

  const String points = polygonMatch[1].str();

  const std::regex pointRegex(R"regex(([-+0-9.eE]+),([-+0-9.eE]+))regex");

  auto begin = std::sregex_iterator(points.begin(), points.end(), pointRegex);

  const auto end = std::sregex_iterator();

  if (begin == end)
    return bounds;

  bounds.minX = std::numeric_limits<double>::max();

  bounds.minY = std::numeric_limits<double>::max();

  bounds.maxX = std::numeric_limits<double>::lowest();

  bounds.maxY = std::numeric_limits<double>::lowest();

  for (auto iterator = begin; iterator != end; ++iterator) {

    const double x = std::stod((*iterator)[1].str());

    const double y = std::stod((*iterator)[2].str());

    bounds.minX = std::min(bounds.minX, x);
    bounds.minY = std::min(bounds.minY, y);
    bounds.maxX = std::max(bounds.maxX, x);
    bounds.maxY = std::max(bounds.maxY, y);
  }

  bounds.valid = bounds.maxX > bounds.minX && bounds.maxY > bounds.minY;

  return bounds;
}

struct SymbolDimensions {
  double width = 0.40;
  double height = 0.40;
};

String symbolStem(TopologicalPowerComp &component) {

  const String type = component.type();

  // Grid-forming / grid-following converter symbols are intentionally
  // orientation-independent. The resource files are GFM.svg and GFL.svg.
  if (contains(type, "GFM"))
    return "GFM";

  if (contains(type, "GFL"))
    return "GFL";

  if (contains(type, "SynchronGenerator"))
    return "generator";

  if (contains(type, "Transformer"))
    return "transformer";

  if (contains(type, "PiLine"))
    return "piline";

  if (contains(type, "Capacitor"))
    return "capacitor";

  if (contains(type, "Inductor") || contains(type, "Reactor"))
    return "inductor";

  if (contains(type, "Resistor"))
    return "resistor";

  if (contains(type, "Switch"))
    return "switch";

  return {};
}

std::filesystem::path
firstExisting(const std::vector<std::filesystem::path> &paths) {

  for (const auto &path : paths) {
    if (!path.empty() && std::filesystem::exists(path))
      return path;
  }

  return {};
}

std::filesystem::path componentSymbol(const std::filesystem::path &directory,
                                      TopologicalPowerComp &component,
                                      const SymbolOrientation orientation) {

  const String stem = symbolStem(component);

  if (stem.empty())
    return {};

  // GFM.svg and GFL.svg are intentionally orientation-independent.
  // Keep the topology-dependent terminal routing/orientation, but use the same
  // visual asset for horizontal and vertical placement.
  if (stem == "GFM" || stem == "GFL") {

    return firstExisting(
        {directory / (stem + ".svg"), directory / (lowerCopy(stem) + ".svg")});
  }

  const String suffix =
      orientation == SymbolOrientation::Horizontal ? "_h.svg" : "_v.svg";

  String capitalized = stem;

  if (!capitalized.empty()) {
    capitalized[0] = static_cast<char>(
        std::toupper(static_cast<unsigned char>(capitalized[0])));
  }

  std::vector<std::filesystem::path> candidates{
      directory / (stem + suffix), directory / (capitalized + suffix)};

  if (stem == "piline") {
    candidates.push_back(directory /
                         (orientation == SymbolOrientation::Horizontal
                              ? "PiLine_h.svg"
                              : "PiLine_v.svg"));
  }

  if (orientation == SymbolOrientation::Horizontal) {

    candidates.push_back(
        directory / (stem == "piline" ? "PiLine.svg" : capitalized + ".svg"));

    // Explicit legacy fallback for a newly added switch symbol.
    if (stem == "switch") {
      candidates.push_back(directory / "Switch.svg");
      candidates.push_back(directory / "switch.svg");
    }
  }

  return firstExisting(candidates);
}

std::filesystem::path busbarSymbol(const std::filesystem::path &directory,
                                   const BusbarOrientation orientation) {

  if (orientation == BusbarOrientation::Vertical) {

    return firstExisting({directory / "busbar_v.svg",
                          directory / "BusBar_v.svg",
                          directory / "Busbar_v.svg", directory / "BusBar.svg",
                          directory / "busbar.svg"});
  }

  return firstExisting({directory / "busbar_h.svg", directory / "BusBar_h.svg",
                        directory / "Busbar_h.svg", directory / "BusBar.svg",
                        directory / "busbar.svg"});
}

std::filesystem::path groundSymbol(const std::filesystem::path &directory,
                                   const String &side) {

  return directory / ("ground_" + lowerCopy(side) + ".svg");
}

SymbolDimensions baseDimensions(const String &type) {

  SymbolDimensions result{RenderStyle::componentWidth,
                          RenderStyle::componentHeight};

  if (contains(type, "GFM") || contains(type, "GFL"))
    return {0.8, 0.8};

  if (contains(type, "SynchronGenerator"))
    return {0.40, 0.40};

  if (contains(type, "Transformer"))
    return {0.50, 0.36};

  if (contains(type, "PiLine"))
    return {0.58, 0.22};

  if (contains(type, "Resistor"))
    return {0.46, 0.22};

  if (contains(type, "Switch"))
    return {0.46, 0.46};

  if (contains(type, "Inductor") || contains(type, "Reactor"))
    return {0.46, 0.22};

  if (contains(type, "Capacitor"))
    return {0.38, 0.28};

  return result;
}

SymbolDimensions orientedDimensions(const String &type,
                                    const SymbolOrientation orientation) {

  auto dimensions = baseDimensions(type);

  if (orientation == SymbolOrientation::Vertical) {

    std::swap(dimensions.width, dimensions.height);
  }

  return dimensions;
}

SvgBounds resizedBounds(const SvgBounds &original,
                        const SymbolDimensions &dimensions) {

  constexpr double pointsPerInch = 72.0;

  const double width = dimensions.width * pointsPerInch;

  const double height = dimensions.height * pointsPerInch;

  SvgBounds result;

  result.minX = centerX(original) - 0.5 * width;

  result.maxX = centerX(original) + 0.5 * width;

  result.minY = centerY(original) - 0.5 * height;

  result.maxY = centerY(original) + 0.5 * height;

  result.valid = true;
  return result;
}

Direction normalize(Direction direction) {
  const double magnitude =
      std::sqrt(direction.x * direction.x + direction.y * direction.y);

  if (magnitude <= 1.0e-12)
    return {1.0, 0.0};

  direction.x /= magnitude;
  direction.y /= magnitude;
  return direction;
}

Direction cardinal(Direction direction) {
  direction = normalize(direction);

  if (std::abs(direction.x) >= std::abs(direction.y)) {

    return {direction.x >= 0.0 ? 1.0 : -1.0, 0.0};
  }

  return {0.0, direction.y >= 0.0 ? 1.0 : -1.0};
}

// -----------------------------------------------------------------------------
// SVG symbol embedding
// -----------------------------------------------------------------------------

String createNestedSvg(const std::filesystem::path &symbolPath,
                       const SvgBounds &bounds, const Bool stretch = false) {

  String source = readTextFile(symbolPath);

  forceMonochromeSvg(source);

  const auto svgStart = source.find("<svg");

  const auto openingEnd = source.find('>', svgStart);

  const auto svgEnd = source.rfind("</svg>");

  if (svgStart == String::npos || openingEnd == String::npos ||
      svgEnd == String::npos || svgEnd <= openingEnd) {

    throw std::runtime_error("Malformed SVG symbol: " + symbolPath.string());
  }

  String openingTag = source.substr(svgStart, openingEnd - svgStart);

  openingTag = std::regex_replace(
      openingTag,
      std::regex(
          R"regex(\s+(width|height|x|y|preserveAspectRatio)\s*=\s*("[^"]*"|'[^']*'))regex",
          std::regex::icase),
      "");

  std::ostringstream nested;

  nested << openingTag << " x=\"" << formatReal(bounds.minX) << "\""
         << " y=\"" << formatReal(bounds.minY) << "\""
         << " width=\"" << formatReal(bounds.width()) << "\""
         << " height=\"" << formatReal(bounds.height()) << "\""
         << " preserveAspectRatio=\"" << (stretch ? "none" : "xMidYMid meet")
         << "\">" << source.substr(openingEnd + 1, svgEnd - openingEnd - 1)
         << "</svg>";

  return nested.str();
}

Bool injectSymbolIntoNode(String &svg, const String &graphNodeName,
                          const std::filesystem::path &symbolPath,
                          const Bool stretch = false) {

  if (!std::filesystem::exists(symbolPath))
    return false;

  const String title = "<title>" + graphNodeName + "</title>";

  const auto titlePosition = svg.find(title);

  if (titlePosition == String::npos)
    return false;

  SvgBounds bounds;

  if (!extractNodeBounds(svg, graphNodeName, bounds))
    return false;

  const String nestedSymbol =
      "\n" + createNestedSvg(symbolPath, bounds, stretch) + "\n";

  // Graphviz represents SVG tooltips using an <a xlink:title="..."> wrapper.
  //
  // Example:
  //
  //   <g class="node">
  //     <title>component_R1</title>
  //     <g id="a_node1">
  //       <a xlink:title="Attributes: ...">
  //         ...
  //       </a>
  //     </g>
  //   </g>
  //
  // The custom electrical SVG must be inserted INSIDE this anchor. Otherwise
  // the transparent Graphviz placeholder has a tooltip but the visible custom
  // symbol itself can intercept the pointer and the tooltip appears unreliable.
  //
  // If no tooltip anchor exists, preserve the previous insertion behaviour.
  const auto searchStart = titlePosition + title.size();

  auto nodeEnd = svg.find("\n<!--", searchStart);

  if (nodeEnd == String::npos) {

    nodeEnd = svg.size();
  }

  const auto anchorStart = svg.find("<a ", searchStart);

  if (anchorStart != String::npos && anchorStart < nodeEnd) {

    const auto anchorOpeningEnd = svg.find('>', anchorStart);

    if (anchorOpeningEnd != String::npos && anchorOpeningEnd < nodeEnd) {

      svg.insert(anchorOpeningEnd + 1, nestedSymbol);

      return true;
    }
  }

  svg.insert(searchStart, nestedSymbol);

  return true;
}

Bool injectFreeSymbol(String &svg, const std::filesystem::path &symbolPath,
                      const SvgBounds &bounds, const String &className) {

  if (!std::filesystem::exists(symbolPath))
    return false;

  const auto graphEnd = svg.rfind("</g>");

  if (graphEnd == String::npos)
    return false;

  svg.insert(graphEnd, "\n<g class=\"" + className + "\">\n" +
                           createNestedSvg(symbolPath, bounds) + "\n</g>\n");

  return true;
}

// -----------------------------------------------------------------------------
// Electrical routing
// -----------------------------------------------------------------------------

enum class BroadSide { Negative, Positive };

struct Connection {
  String bus;
  String component;

  SvgBounds busBounds;
  SvgBounds componentBounds;

  BusbarOrientation busOrientation = BusbarOrientation::Vertical;

  SymbolOrientation componentOrientation = SymbolOrientation::Horizontal;

  BroadSide side = BroadSide::Negative;

  Bool singleBusComponent = false;
  Bool groundedLeaf = false;
};

String edgeTitle(const String &first, const String &second) {

  return "<title>" + first + "&#45;&gt;" + second + "</title>";
}

Bool replaceEdgePath(String &svg, const String &first, const String &second,
                     const String &pathData) {

  String title = edgeTitle(first, second);

  auto titlePosition = svg.find(title);

  if (titlePosition == String::npos) {
    title = edgeTitle(second, first);
    titlePosition = svg.find(title);
  }

  if (titlePosition == String::npos)
    return false;

  const auto groupEnd = svg.find("</g>", titlePosition);

  if (groupEnd == String::npos)
    return false;

  String group = svg.substr(titlePosition, groupEnd + 4 - titlePosition);

  const std::regex pathRegex(R"regex(<path([^>]*)\bd="([^"]+)"([^>]*)/>)regex");

  std::smatch match;

  if (!std::regex_search(group, match, pathRegex))
    return false;

  const String replacement = "<path" + match[1].str() + "d=\"" + pathData +
                             "\"" + match[3].str() + "/>";

  group.replace(static_cast<std::size_t>(match.position()),
                static_cast<std::size_t>(match.length()), replacement);

  svg.replace(titlePosition, groupEnd + 4 - titlePosition, group);

  return true;
}

BroadSide broadSideFor(const SvgBounds &bus, const SvgBounds &component,
                       const BusbarOrientation orientation) {

  if (orientation == BusbarOrientation::Vertical) {

    return centerX(component) < centerX(bus) ? BroadSide::Negative
                                             : BroadSide::Positive;
  }

  return centerY(component) < centerY(bus) ? BroadSide::Negative
                                           : BroadSide::Positive;
}

double longCoordinate(const SvgBounds &bounds,
                      const BusbarOrientation orientation) {

  return orientation == BusbarOrientation::Vertical ? centerY(bounds)
                                                    : centerX(bounds);
}

double busLongMin(const SvgBounds &bounds,
                  const BusbarOrientation orientation) {

  return orientation == BusbarOrientation::Vertical ? bounds.minY : bounds.minX;
}

double busLongLength(const SvgBounds &bounds,
                     const BusbarOrientation orientation) {

  return orientation == BusbarOrientation::Vertical ? bounds.height()
                                                    : bounds.width();
}

struct RoutePoint {
  double x = 0.0;
  double y = 0.0;
};

struct RouteSegment {
  RoutePoint first;
  RoutePoint second;
};

struct RouteObstacle {
  String name;
  SvgBounds bounds;
};

struct RoutePath {
  String bus;
  String component;
  std::vector<RoutePoint> points;
  double score = std::numeric_limits<double>::max();
};

struct RouteBridge {
  std::size_t segmentIndex = 0;
  double coordinate = 0.0;
};

Bool samePoint(const RoutePoint &first, const RoutePoint &second,
               const double tolerance = 1.0e-6) {

  return std::abs(first.x - second.x) <= tolerance &&
         std::abs(first.y - second.y) <= tolerance;
}

void appendRoutePoint(std::vector<RoutePoint> &points,
                      const RoutePoint &point) {

  if (!points.empty() && samePoint(points.back(), point)) {

    return;
  }

  points.push_back(point);
}

Bool routeSegmentHorizontal(const RouteSegment &segment) {

  return std::abs(segment.first.y - segment.second.y) < 1.0e-7;
}

Bool routeSegmentVertical(const RouteSegment &segment) {

  return std::abs(segment.first.x - segment.second.x) < 1.0e-7;
}

std::vector<RoutePoint>
simplifyRoutePoints(const std::vector<RoutePoint> &input) {

  std::vector<RoutePoint> points;

  for (const auto &point : input) {

    appendRoutePoint(points, point);
  }

  if (points.size() < 3) {

    return points;
  }

  Bool changed = true;

  while (changed && points.size() >= 3) {

    changed = false;

    std::vector<RoutePoint> simplified;

    simplified.push_back(points.front());

    for (std::size_t index = 1; index + 1 < points.size(); ++index) {

      const RoutePoint &previous = simplified.back();

      const RoutePoint &current = points[index];

      const RoutePoint &next = points[index + 1];

      const Bool horizontal = std::abs(previous.y - current.y) < 1.0e-7 &&
                              std::abs(current.y - next.y) < 1.0e-7;

      const Bool vertical = std::abs(previous.x - current.x) < 1.0e-7 &&
                            std::abs(current.x - next.x) < 1.0e-7;

      if (horizontal || vertical) {

        changed = true;

        continue;
      }

      simplified.push_back(current);
    }

    simplified.push_back(points.back());

    points = std::move(simplified);
  }

  return points;
}

std::vector<RouteSegment> routeSegments(const std::vector<RoutePoint> &points) {

  std::vector<RouteSegment> segments;

  if (points.size() < 2) {

    return segments;
  }

  segments.reserve(points.size() - 1);

  for (std::size_t index = 1; index < points.size(); ++index) {

    segments.push_back({points[index - 1], points[index]});
  }

  return segments;
}

double routeLength(const std::vector<RoutePoint> &points) {

  double length = 0.0;

  for (const auto &segment : routeSegments(points)) {

    length += std::abs(segment.second.x - segment.first.x) +
              std::abs(segment.second.y - segment.first.y);
  }

  return length;
}

Bool routeSegmentIntersectsBounds(const RouteSegment &segment, SvgBounds bounds,
                                  const double clearance) {

  bounds.minX -= clearance;

  bounds.maxX += clearance;

  bounds.minY -= clearance;

  bounds.maxY += clearance;

  const double minX = std::min(segment.first.x, segment.second.x);

  const double maxX = std::max(segment.first.x, segment.second.x);

  const double minY = std::min(segment.first.y, segment.second.y);

  const double maxY = std::max(segment.first.y, segment.second.y);

  if (routeSegmentHorizontal(segment)) {

    return segment.first.y >= bounds.minY && segment.first.y <= bounds.maxY &&
           maxX >= bounds.minX && minX <= bounds.maxX;
  }

  if (routeSegmentVertical(segment)) {

    return segment.first.x >= bounds.minX && segment.first.x <= bounds.maxX &&
           maxY >= bounds.minY && minY <= bounds.maxY;
  }

  return maxX >= bounds.minX && minX <= bounds.maxX && maxY >= bounds.minY &&
         minY <= bounds.maxY;
}

Bool routeProperCrossing(const RouteSegment &first, const RouteSegment &second,
                         RoutePoint *crossing = nullptr) {

  constexpr double tolerance = 1.0e-7;

  const Bool firstHorizontal = routeSegmentHorizontal(first);

  const Bool firstVertical = routeSegmentVertical(first);

  const Bool secondHorizontal = routeSegmentHorizontal(second);

  const Bool secondVertical = routeSegmentVertical(second);

  if (!((firstHorizontal && secondVertical) ||
        (firstVertical && secondHorizontal))) {

    return false;
  }

  const RouteSegment &horizontal = firstHorizontal ? first : second;

  const RouteSegment &vertical = firstVertical ? first : second;

  const double horizontalMinX =
      std::min(horizontal.first.x, horizontal.second.x);

  const double horizontalMaxX =
      std::max(horizontal.first.x, horizontal.second.x);

  const double verticalMinY = std::min(vertical.first.y, vertical.second.y);

  const double verticalMaxY = std::max(vertical.first.y, vertical.second.y);

  const double x = vertical.first.x;

  const double y = horizontal.first.y;

  if (x <= horizontalMinX + tolerance || x >= horizontalMaxX - tolerance ||
      y <= verticalMinY + tolerance || y >= verticalMaxY - tolerance) {

    return false;
  }

  if (crossing) {

    crossing->x = x;

    crossing->y = y;
  }

  return true;
}

Bool routeCollinearOverlap(const RouteSegment &first,
                           const RouteSegment &second) {

  constexpr double tolerance = 1.0e-7;

  if (routeSegmentHorizontal(first) && routeSegmentHorizontal(second) &&
      std::abs(first.first.y - second.first.y) <= tolerance) {

    const double overlap = std::min(std::max(first.first.x, first.second.x),
                                    std::max(second.first.x, second.second.x)) -
                           std::max(std::min(first.first.x, first.second.x),
                                    std::min(second.first.x, second.second.x));

    return overlap > 0.5;
  }

  if (routeSegmentVertical(first) && routeSegmentVertical(second) &&
      std::abs(first.first.x - second.first.x) <= tolerance) {

    const double overlap = std::min(std::max(first.first.y, first.second.y),
                                    std::max(second.first.y, second.second.y)) -
                           std::max(std::min(first.first.y, first.second.y),
                                    std::min(second.first.y, second.second.y));

    return overlap > 0.5;
  }

  return false;
}

std::vector<RouteObstacle>
routeObstacles(const std::vector<Connection> &connections) {

  std::vector<RouteObstacle> obstacles;

  std::unordered_set<String> seen;

  for (const auto &connection : connections) {

    if (seen.insert(connection.bus).second) {

      obstacles.push_back({connection.bus, connection.busBounds});
    }

    if (seen.insert(connection.component).second) {

      obstacles.push_back({connection.component, connection.componentBounds});
    }
  }

  return obstacles;
}

struct RouteEndpoints {
  RoutePoint bus;
  RoutePoint escape;
  RoutePoint approach;
  RoutePoint terminal;

  Bool verticalBus = false;

  double outwardSign = 1.0;
};

RouteEndpoints routeEndpoints(const Connection &connection, const double port) {

  constexpr double busbarEscape = 8.0;

  constexpr double terminalApproach = 10.0;

  RouteEndpoints result;

  result.verticalBus = connection.busOrientation == BusbarOrientation::Vertical;

  result.outwardSign = connection.side == BroadSide::Negative ? -1.0 : 1.0;

  if (result.verticalBus) {

    result.bus = {connection.side == BroadSide::Negative
                      ? connection.busBounds.minX
                      : connection.busBounds.maxX,
                  port};

    result.escape = {result.bus.x + result.outwardSign * busbarEscape,
                     result.bus.y};

  } else {

    result.bus = {port, connection.side == BroadSide::Negative
                            ? connection.busBounds.minY
                            : connection.busBounds.maxY};

    result.escape = {result.bus.x,
                     result.bus.y + result.outwardSign * busbarEscape};
  }

  if (connection.componentOrientation == SymbolOrientation::Horizontal) {

    const Bool componentOnRight =
        centerX(connection.componentBounds) >= centerX(connection.busBounds);

    result.terminal = {componentOnRight ? connection.componentBounds.minX
                                        : connection.componentBounds.maxX,
                       centerY(connection.componentBounds)};

    result.approach = {componentOnRight ? result.terminal.x - terminalApproach
                                        : result.terminal.x + terminalApproach,
                       result.terminal.y};

  } else {

    const Bool componentBelow =
        centerY(connection.componentBounds) >= centerY(connection.busBounds);

    result.terminal = {centerX(connection.componentBounds),
                       componentBelow ? connection.componentBounds.minY
                                      : connection.componentBounds.maxY};

    result.approach = {result.terminal.x,
                       componentBelow ? result.terminal.y - terminalApproach
                                      : result.terminal.y + terminalApproach};
  }

  return result;
}

void addRouteCandidate(std::vector<RoutePath> &candidates,
                       const Connection &connection,
                       const RouteEndpoints &endpoints,
                       const std::vector<RoutePoint> &middle) {

  RoutePath candidate;

  candidate.bus = connection.bus;

  candidate.component = connection.component;

  appendRoutePoint(candidate.points, endpoints.bus);

  appendRoutePoint(candidate.points, endpoints.escape);

  for (const auto &point : middle) {

    appendRoutePoint(candidate.points, point);
  }

  appendRoutePoint(candidate.points, endpoints.approach);

  appendRoutePoint(candidate.points, endpoints.terminal);

  candidate.points = simplifyRoutePoints(candidate.points);

  candidates.push_back(std::move(candidate));
}

std::vector<RoutePath> routeCandidates(const Connection &connection,
                                       const double port) {

  constexpr double laneSpacing = 18.0;

  constexpr UInt laneCount = 14;

  const RouteEndpoints endpoints = routeEndpoints(connection, port);

  std::vector<RoutePath> candidates;

  // Direct middle connection where possible.
  if (std::abs(endpoints.escape.x - endpoints.approach.x) < 1.0e-7 ||
      std::abs(endpoints.escape.y - endpoints.approach.y) < 1.0e-7) {

    addRouteCandidate(candidates, connection, endpoints, {});
  }

  // Standard Manhattan candidates.
  addRouteCandidate(candidates, connection, endpoints,
                    {{endpoints.approach.x, endpoints.escape.y}});

  addRouteCandidate(candidates, connection, endpoints,
                    {{endpoints.escape.x, endpoints.approach.y}});

  // Additional routing lanes. These are what make the checker useful on dense
  // systems: if the short route collides with an object/edge, the connection is
  // allowed to use a more distant lane and therefore enlarge the final SVG.
  for (UInt lane = 1; lane <= laneCount; ++lane) {

    const double offset = static_cast<double>(lane) * laneSpacing;

    const double upperY =
        std::min(endpoints.escape.y, endpoints.approach.y) - offset;

    const double lowerY =
        std::max(endpoints.escape.y, endpoints.approach.y) + offset;

    addRouteCandidate(
        candidates, connection, endpoints,
        {{endpoints.escape.x, upperY}, {endpoints.approach.x, upperY}});

    addRouteCandidate(
        candidates, connection, endpoints,
        {{endpoints.escape.x, lowerY}, {endpoints.approach.x, lowerY}});

    const double leftX =
        std::min(endpoints.escape.x, endpoints.approach.x) - offset;

    const double rightX =
        std::max(endpoints.escape.x, endpoints.approach.x) + offset;

    addRouteCandidate(
        candidates, connection, endpoints,
        {{leftX, endpoints.escape.y}, {leftX, endpoints.approach.y}});

    addRouteCandidate(
        candidates, connection, endpoints,
        {{rightX, endpoints.escape.y}, {rightX, endpoints.approach.y}});
  }

  return candidates;
}

double scoreRoute(const RoutePath &candidate,
                  const std::vector<RouteObstacle> &obstacles,
                  const std::vector<RoutePath> &existingRoutes) {

  constexpr double objectClearance = 4.0;

  constexpr double objectPenalty = 1.0e8;

  constexpr double collinearPenalty = 1.0e7;

  constexpr double crossingPenalty = 1.0e5;

  double score = routeLength(candidate.points);

  if (candidate.points.size() > 2) {

    score += 3.0 * static_cast<double>(candidate.points.size() - 2);
  }

  const auto candidateSegments = routeSegments(candidate.points);

  for (std::size_t segmentIndex = 0; segmentIndex < candidateSegments.size();
       ++segmentIndex) {

    const auto &segment = candidateSegments[segmentIndex];

    for (const auto &obstacle : obstacles) {

      if (obstacle.name == candidate.bus && segmentIndex == 0)
        continue;

      if (obstacle.name == candidate.component &&
          segmentIndex + 1 == candidateSegments.size())
        continue;

      if (routeSegmentIntersectsBounds(segment, obstacle.bounds,
                                       objectClearance)) {

        score += objectPenalty;
      }
    }

    for (const auto &existingRoute : existingRoutes) {

      for (const auto &existingSegment : routeSegments(existingRoute.points)) {

        if (routeCollinearOverlap(segment, existingSegment)) {

          score += collinearPenalty;

        } else if (routeProperCrossing(segment, existingSegment)) {

          score += crossingPenalty;
        }
      }
    }
  }

  return score;
}

RoutePath chooseBestRoute(const Connection &connection, const double port,
                          const std::vector<RouteObstacle> &obstacles,
                          const std::vector<RoutePath> &existingRoutes) {

  auto candidates = routeCandidates(connection, port);

  RoutePath best;

  for (auto &candidate : candidates) {

    candidate.score = scoreRoute(candidate, obstacles, existingRoutes);

    if (candidate.score < best.score) {

      best = candidate;
    }
  }

  return best;
}

String routePathWithBridges(const RoutePath &route,
                            std::vector<RouteBridge> bridges) {

  constexpr double bridgeRadius = 4.0;

  constexpr double bridgeHeight = 5.0;

  if (route.points.empty()) {
    return {};
  }

  std::ostringstream path;

  path << "M" << formatReal(route.points.front().x) << ","
       << formatReal(route.points.front().y);

  const auto segments = routeSegments(route.points);

  for (std::size_t segmentIndex = 0; segmentIndex < segments.size();
       ++segmentIndex) {

    const auto &segment = segments[segmentIndex];

    std::vector<RouteBridge> segmentBridges;

    for (const auto &bridge : bridges) {

      if (bridge.segmentIndex == segmentIndex) {

        segmentBridges.push_back(bridge);
      }
    }

    if (segmentBridges.empty() ||
        (!routeSegmentHorizontal(segment) && !routeSegmentVertical(segment))) {

      path << " L" << formatReal(segment.second.x) << ","
           << formatReal(segment.second.y);

      continue;
    }

    const Bool horizontal = routeSegmentHorizontal(segment);

    const double direction =
        horizontal ? (segment.second.x >= segment.first.x ? 1.0 : -1.0)
                   : (segment.second.y >= segment.first.y ? 1.0 : -1.0);

    std::sort(segmentBridges.begin(), segmentBridges.end(),
              [&](const RouteBridge &first, const RouteBridge &second) {
                return direction > 0.0 ? first.coordinate < second.coordinate
                                       : first.coordinate > second.coordinate;
              });

    for (const auto &bridge : segmentBridges) {

      if (horizontal) {

        const double preX = bridge.coordinate - direction * bridgeRadius;

        const double postX = bridge.coordinate + direction * bridgeRadius;

        if ((direction > 0.0 &&
             (preX <= segment.first.x || postX >= segment.second.x)) ||
            (direction < 0.0 &&
             (preX >= segment.first.x || postX <= segment.second.x))) {

          continue;
        }

        path << " L" << formatReal(preX) << ","
             << formatReal(segment.first.y)

             // Small bridge/hump over the crossing line.
             << " Q" << formatReal(bridge.coordinate) << ","
             << formatReal(segment.first.y - bridgeHeight) << " "
             << formatReal(postX) << "," << formatReal(segment.first.y);

      } else {

        const double preY = bridge.coordinate - direction * bridgeRadius;

        const double postY = bridge.coordinate + direction * bridgeRadius;

        if ((direction > 0.0 &&
             (preY <= segment.first.y || postY >= segment.second.y)) ||
            (direction < 0.0 &&
             (preY >= segment.first.y || postY <= segment.second.y))) {

          continue;
        }

        path << " L" << formatReal(segment.first.x) << "," << formatReal(preY)

             << " Q" << formatReal(segment.first.x + bridgeHeight) << ","
             << formatReal(bridge.coordinate) << " "
             << formatReal(segment.first.x) << "," << formatReal(postY);
      }
    }

    path << " L" << formatReal(segment.second.x) << ","
         << formatReal(segment.second.y);
  }

  return path.str();
}

struct OrthogonalQueueItem {
  double priority = 0.0;

  std::size_t state = 0;

  Bool operator>(const OrthogonalQueueItem &other) const {

    return priority > other.priority;
  }
};

void addRoutingCoordinate(std::vector<double> &coordinates,
                          const double value) {

  coordinates.push_back(value);
}

void normalizeRoutingCoordinates(std::vector<double> &coordinates) {

  std::sort(coordinates.begin(), coordinates.end());

  std::vector<double> unique;

  constexpr double tolerance = 0.25;

  for (const double value : coordinates) {

    if (unique.empty() || std::abs(value - unique.back()) > tolerance) {

      unique.push_back(value);
    }
  }

  coordinates = std::move(unique);
}

std::size_t routingCoordinateIndex(const std::vector<double> &coordinates,
                                   const double value) {

  auto iterator =
      std::lower_bound(coordinates.begin(), coordinates.end(), value);

  if (iterator == coordinates.end()) {

    return coordinates.size() - 1;
  }

  if (iterator != coordinates.begin()) {

    const auto previous = std::prev(iterator);

    if (std::abs(*previous - value) < std::abs(*iterator - value)) {

      iterator = previous;
    }
  }

  return static_cast<std::size_t>(std::distance(coordinates.begin(), iterator));
}

struct RoutingWorkspace {
  std::vector<double> xCoordinates;
  std::vector<double> yCoordinates;
  std::vector<unsigned char> horizontalBlocked;
  std::vector<unsigned char> verticalBlocked;
  std::size_t xCount = 0;
  std::size_t yCount = 0;
};

struct WorkspaceEndpointPair {
  RoutePoint escape;
  RoutePoint approach;
};

std::size_t horizontalWorkspaceEdgeIndex(const RoutingWorkspace &workspace,
                                         const std::size_t x,
                                         const std::size_t y) {
  return y * (workspace.xCount - 1) + x;
}

std::size_t verticalWorkspaceEdgeIndex(const RoutingWorkspace &workspace,
                                       const std::size_t x,
                                       const std::size_t y) {
  return x * (workspace.yCount - 1) + y;
}

Bool workspaceEdgeBlocked(const RoutingWorkspace &workspace,
                          const std::size_t x0, const std::size_t y0,
                          const std::size_t x1, const std::size_t y1) {

  if (y0 == y1 && x0 + 1 == x1)
    return workspace.horizontalBlocked[horizontalWorkspaceEdgeIndex(
               workspace, x0, y0)] != 0;
  if (y0 == y1 && x1 + 1 == x0)
    return workspace.horizontalBlocked[horizontalWorkspaceEdgeIndex(
               workspace, x1, y0)] != 0;
  if (x0 == x1 && y0 + 1 == y1)
    return workspace.verticalBlocked[verticalWorkspaceEdgeIndex(workspace, x0,
                                                                y0)] != 0;
  if (x0 == x1 && y1 + 1 == y0)
    return workspace.verticalBlocked[verticalWorkspaceEdgeIndex(workspace, x0,
                                                                y1)] != 0;
  return true;
}

RoutingWorkspace
buildRoutingWorkspace(const std::vector<RouteObstacle> &obstacles,
                      const std::vector<WorkspaceEndpointPair> &endpointPairs) {

  constexpr double obstacleTrackClearance = 9.0;
  constexpr double objectClearance = 4.0;
  RoutingWorkspace workspace;

  double minimumX = std::numeric_limits<double>::max();
  double maximumX = std::numeric_limits<double>::lowest();
  double minimumY = std::numeric_limits<double>::max();
  double maximumY = std::numeric_limits<double>::lowest();

  auto includePoint = [&](const RoutePoint &point) {
    addRoutingCoordinate(workspace.xCoordinates, point.x);
    addRoutingCoordinate(workspace.yCoordinates, point.y);
    minimumX = std::min(minimumX, point.x);
    maximumX = std::max(maximumX, point.x);
    minimumY = std::min(minimumY, point.y);
    maximumY = std::max(maximumY, point.y);
  };

  for (const auto &pair : endpointPairs) {
    includePoint(pair.escape);
    includePoint(pair.approach);
  }

  for (const auto &obstacle : obstacles) {
    addRoutingCoordinate(workspace.xCoordinates,
                         obstacle.bounds.minX - obstacleTrackClearance);
    addRoutingCoordinate(workspace.xCoordinates,
                         obstacle.bounds.maxX + obstacleTrackClearance);
    addRoutingCoordinate(workspace.yCoordinates,
                         obstacle.bounds.minY - obstacleTrackClearance);
    addRoutingCoordinate(workspace.yCoordinates,
                         obstacle.bounds.maxY + obstacleTrackClearance);
    minimumX = std::min(minimumX, obstacle.bounds.minX);
    maximumX = std::max(maximumX, obstacle.bounds.maxX);
    minimumY = std::min(minimumY, obstacle.bounds.minY);
    maximumY = std::max(maximumY, obstacle.bounds.maxY);
  }

  if (minimumX > maximumX) {
    minimumX = -100.0;
    maximumX = 100.0;
    minimumY = -100.0;
    maximumY = 100.0;
  }

  // Global outer lanes mean difficult routes can enlarge the drawing without
  // rebuilding the routing grid per connection.
  for (const double margin : {70.0, 150.0, 260.0}) {
    addRoutingCoordinate(workspace.xCoordinates, minimumX - margin);
    addRoutingCoordinate(workspace.xCoordinates, maximumX + margin);
    addRoutingCoordinate(workspace.yCoordinates, minimumY - margin);
    addRoutingCoordinate(workspace.yCoordinates, maximumY + margin);
  }

  normalizeRoutingCoordinates(workspace.xCoordinates);
  normalizeRoutingCoordinates(workspace.yCoordinates);
  workspace.xCount = workspace.xCoordinates.size();
  workspace.yCount = workspace.yCoordinates.size();

  if (workspace.xCount < 2 || workspace.yCount < 2)
    return workspace;

  workspace.horizontalBlocked.assign(workspace.yCount * (workspace.xCount - 1),
                                     0);
  workspace.verticalBlocked.assign(workspace.xCount * (workspace.yCount - 1),
                                   0);

  // Static obstacle tests are paid ONCE for the complete rendering pass.
  for (std::size_t y = 0; y < workspace.yCount; ++y) {
    for (std::size_t x = 0; x + 1 < workspace.xCount; ++x) {
      const RouteSegment segment{
          {workspace.xCoordinates[x], workspace.yCoordinates[y]},
          {workspace.xCoordinates[x + 1], workspace.yCoordinates[y]}};
      for (const auto &obstacle : obstacles) {
        if (routeSegmentIntersectsBounds(segment, obstacle.bounds,
                                         objectClearance)) {
          workspace.horizontalBlocked[horizontalWorkspaceEdgeIndex(workspace, x,
                                                                   y)] = 1;
          break;
        }
      }
    }
  }

  for (std::size_t x = 0; x < workspace.xCount; ++x) {
    for (std::size_t y = 0; y + 1 < workspace.yCount; ++y) {
      const RouteSegment segment{
          {workspace.xCoordinates[x], workspace.yCoordinates[y]},
          {workspace.xCoordinates[x], workspace.yCoordinates[y + 1]}};
      for (const auto &obstacle : obstacles) {
        if (routeSegmentIntersectsBounds(segment, obstacle.bounds,
                                         objectClearance)) {
          workspace
              .verticalBlocked[verticalWorkspaceEdgeIndex(workspace, x, y)] = 1;
          break;
        }
      }
    }
  }

  return workspace;
}

struct RoutingOccupancy {
  std::vector<unsigned char> horizontalUsed;
  std::vector<unsigned char> verticalUsed;
  std::vector<unsigned char> nodeHorizontalUsed;
  std::vector<unsigned char> nodeVerticalUsed;
};

RoutingOccupancy makeRoutingOccupancy(const RoutingWorkspace &workspace) {

  RoutingOccupancy occupancy;
  occupancy.horizontalUsed.assign(workspace.horizontalBlocked.size(), 0);
  occupancy.verticalUsed.assign(workspace.verticalBlocked.size(), 0);
  occupancy.nodeHorizontalUsed.assign(workspace.xCount * workspace.yCount, 0);
  occupancy.nodeVerticalUsed.assign(workspace.xCount * workspace.yCount, 0);
  return occupancy;
}

std::size_t workspaceNodeIndex(const RoutingWorkspace &workspace,
                               const std::size_t x, const std::size_t y) {
  return x * workspace.yCount + y;
}

Bool exactRoutingCoordinate(const std::vector<double> &coordinates,
                            const double value, std::size_t &index) {

  if (coordinates.empty())
    return false;
  index = routingCoordinateIndex(coordinates, value);
  return std::abs(coordinates[index] - value) <= 1.0e-4;
}

void markRouteInOccupancy(const RoutingWorkspace &workspace,
                          RoutingOccupancy &occupancy, const RoutePath &route) {

  for (const auto &segment : routeSegments(route.points)) {
    std::size_t x0 = 0, y0 = 0, x1 = 0, y1 = 0;
    if (!exactRoutingCoordinate(workspace.xCoordinates, segment.first.x, x0) ||
        !exactRoutingCoordinate(workspace.yCoordinates, segment.first.y, y0) ||
        !exactRoutingCoordinate(workspace.xCoordinates, segment.second.x, x1) ||
        !exactRoutingCoordinate(workspace.yCoordinates, segment.second.y, y1))
      continue;

    if (y0 == y1 && x0 != x1) {
      const std::size_t begin = std::min(x0, x1);
      const std::size_t end = std::max(x0, x1);
      for (std::size_t x = begin; x < end; ++x) {
        occupancy
            .horizontalUsed[horizontalWorkspaceEdgeIndex(workspace, x, y0)] = 1;
        occupancy.nodeHorizontalUsed[workspaceNodeIndex(workspace, x, y0)] = 1;
        occupancy.nodeHorizontalUsed[workspaceNodeIndex(workspace, x + 1, y0)] =
            1;
      }
    } else if (x0 == x1 && y0 != y1) {
      const std::size_t begin = std::min(y0, y1);
      const std::size_t end = std::max(y0, y1);
      for (std::size_t y = begin; y < end; ++y) {
        occupancy.verticalUsed[verticalWorkspaceEdgeIndex(workspace, x0, y)] =
            1;
        occupancy.nodeVerticalUsed[workspaceNodeIndex(workspace, x0, y)] = 1;
        occupancy.nodeVerticalUsed[workspaceNodeIndex(workspace, x0, y + 1)] =
            1;
      }
    }
  }
}

double workspaceWireCost(const RoutingWorkspace &workspace,
                         const RoutingOccupancy &occupancy,
                         const std::size_t x0, const std::size_t y0,
                         const std::size_t x1, const std::size_t y1) {

  constexpr double crossingCost = 550.0;
  constexpr double collinearCost = 50000.0;

  if (y0 == y1) {
    const std::size_t x = std::min(x0, x1);
    double cost =
        occupancy.horizontalUsed[horizontalWorkspaceEdgeIndex(workspace, x, y0)]
            ? collinearCost
            : 0.0;
    if (occupancy.nodeVerticalUsed[workspaceNodeIndex(workspace, x0, y0)] ||
        occupancy.nodeVerticalUsed[workspaceNodeIndex(workspace, x1, y1)])
      cost += crossingCost;
    return cost;
  }

  const std::size_t y = std::min(y0, y1);
  double cost =
      occupancy.verticalUsed[verticalWorkspaceEdgeIndex(workspace, x0, y)]
          ? collinearCost
          : 0.0;
  if (occupancy.nodeHorizontalUsed[workspaceNodeIndex(workspace, x0, y0)] ||
      occupancy.nodeHorizontalUsed[workspaceNodeIndex(workspace, x1, y1)])
    cost += crossingCost;
  return cost;
}

Bool routeIsCollisionFree(const RoutePath &route, const Connection &connection,
                          const std::vector<RouteObstacle> &obstacles,
                          const std::vector<RouteSegment> &occupiedSegments) {

  constexpr double objectClearance = 4.0;
  const auto segments = routeSegments(route.points);

  for (std::size_t segmentIndex = 0; segmentIndex < segments.size();
       ++segmentIndex) {
    const auto &segment = segments[segmentIndex];

    for (const auto &obstacle : obstacles) {
      // The own busbar is exempt ONLY on the first outward segment.
      // All later segments must treat the whole busbar, including the short
      // ends, as a hard obstacle.
      if (obstacle.name == connection.bus && segmentIndex == 0)
        continue;

      // Only the final terminal approach may enter the own component.
      if (obstacle.name == connection.component &&
          segmentIndex + 1 == segments.size())
        continue;

      if (routeSegmentIntersectsBounds(segment, obstacle.bounds,
                                       objectClearance))
        return false;
    }

    for (const auto &occupied : occupiedSegments) {
      if (routeCollinearOverlap(segment, occupied) ||
          routeProperCrossing(segment, occupied))
        return false;
    }
  }
  return true;
}

Bool routeStartsAtLongBusSide(const Connection &connection,
                              const RoutePath &route, const double port) {

  constexpr double tolerance = 1.0e-4;
  if (route.points.size() < 2)
    return false;

  const RouteEndpoints expected = routeEndpoints(connection, port);
  if (!samePoint(route.points.front(), expected.bus, tolerance))
    return false;

  const RoutePoint &second = route.points[1];

  if (connection.busOrientation == BusbarOrientation::Vertical) {
    // Vertical busbar: ONLY LEFT/RIGHT long sides. First segment MUST leave
    // horizontally and outward. Never connect at top/bottom short ends.
    if (std::abs(second.y - expected.bus.y) > tolerance)
      return false;
    if (expected.outwardSign * (second.x - expected.bus.x) <= tolerance)
      return false;
    return port > connection.busBounds.minY + tolerance &&
           port < connection.busBounds.maxY - tolerance;
  }

  // Horizontal busbar: ONLY TOP/BOTTOM long sides. First segment MUST leave
  // vertically and outward. Never connect at left/right short ends.
  if (std::abs(second.x - expected.bus.x) > tolerance)
    return false;
  if (expected.outwardSign * (second.y - expected.bus.y) <= tolerance)
    return false;
  return port > connection.busBounds.minX + tolerance &&
         port < connection.busBounds.maxX - tolerance;
}

RoutePath
cheapOrthogonalRoute(const Connection &connection, const double port,
                     const std::vector<RouteObstacle> &obstacles,
                     const std::vector<RouteSegment> &occupiedSegments) {

  const RouteEndpoints endpoints = routeEndpoints(connection, port);

  std::vector<RoutePath> candidates;

  // Straight middle section.
  if (std::abs(endpoints.escape.x - endpoints.approach.x) < 1.0e-7 ||
      std::abs(endpoints.escape.y - endpoints.approach.y) < 1.0e-7) {

    addRouteCandidate(candidates, connection, endpoints, {});
  }

  // Horizontal-first Manhattan connection.
  addRouteCandidate(candidates, connection, endpoints,
                    {{endpoints.approach.x, endpoints.escape.y}});

  // Vertical-first Manhattan connection.
  addRouteCandidate(candidates, connection, endpoints,
                    {{endpoints.escape.x, endpoints.approach.y}});

  RoutePath best;

  for (auto &candidate : candidates) {

    if (!routeIsCollisionFree(candidate, connection, obstacles,
                              occupiedSegments)) {

      continue;
    }

    candidate.score = routeLength(candidate.points) +
                      3.0 * static_cast<double>(candidate.points.size());

    if (candidate.score < best.score) {

      best = std::move(candidate);
    }
  }

  return best;
}

RoutePath orthogonalGridRoute(const Connection &connection, const double port,
                              const RoutingWorkspace &workspace,
                              const RoutingOccupancy &occupancy) {

  const RouteEndpoints endpoints = routeEndpoints(connection, port);
  if (workspace.xCount < 2 || workspace.yCount < 2)
    return {};

  constexpr double turnCost = 9.0;
  constexpr std::size_t directionCount = 3;
  const std::size_t xCount = workspace.xCount;
  const std::size_t yCount = workspace.yCount;
  const std::size_t stateCount = xCount * yCount * directionCount;

  const auto encode = [&](std::size_t x, std::size_t y, std::size_t direction) {
    return (x * yCount + y) * directionCount + direction;
  };
  const auto decode = [&](std::size_t state, std::size_t &x, std::size_t &y,
                          std::size_t &direction) {
    direction = state % directionCount;
    const std::size_t pointState = state / directionCount;
    y = pointState % yCount;
    x = pointState / yCount;
  };

  const std::size_t startX =
      routingCoordinateIndex(workspace.xCoordinates, endpoints.escape.x);
  const std::size_t startY =
      routingCoordinateIndex(workspace.yCoordinates, endpoints.escape.y);
  const std::size_t goalX =
      routingCoordinateIndex(workspace.xCoordinates, endpoints.approach.x);
  const std::size_t goalY =
      routingCoordinateIndex(workspace.yCoordinates, endpoints.approach.y);
  const std::size_t startState = encode(startX, startY, 0);

  const double infinity = std::numeric_limits<double>::infinity();
  std::vector<double> distance(stateCount, infinity);
  std::vector<std::size_t> predecessor(stateCount,
                                       std::numeric_limits<std::size_t>::max());
  std::priority_queue<OrthogonalQueueItem, std::vector<OrthogonalQueueItem>,
                      std::greater<OrthogonalQueueItem>>
      queue;

  auto heuristic = [&](std::size_t x, std::size_t y) {
    return std::abs(workspace.xCoordinates[x] - workspace.xCoordinates[goalX]) +
           std::abs(workspace.yCoordinates[y] - workspace.yCoordinates[goalY]);
  };

  distance[startState] = 0.0;
  queue.push({heuristic(startX, startY), startState});
  std::size_t goalState = std::numeric_limits<std::size_t>::max();

  while (!queue.empty()) {
    const OrthogonalQueueItem item = queue.top();
    queue.pop();

    std::size_t x = 0, y = 0, direction = 0;
    decode(item.state, x, y, direction);
    const double expectedPriority = distance[item.state] + heuristic(x, y);
    if (item.priority > expectedPriority + 1.0e-7)
      continue;

    if (x == goalX && y == goalY) {
      goalState = item.state;
      break;
    }

    struct Neighbor {
      std::size_t x = 0;
      std::size_t y = 0;
      std::size_t direction = 0;
    };
    Neighbor candidates[4];
    std::size_t candidateCount = 0;
    if (x > 0)
      candidates[candidateCount++] = {x - 1, y, 1};
    if (x + 1 < xCount)
      candidates[candidateCount++] = {x + 1, y, 1};
    if (y > 0)
      candidates[candidateCount++] = {x, y - 1, 2};
    if (y + 1 < yCount)
      candidates[candidateCount++] = {x, y + 1, 2};

    const RoutePoint currentPoint{workspace.xCoordinates[x],
                                  workspace.yCoordinates[y]};

    for (std::size_t candidateIndex = 0; candidateIndex < candidateCount;
         ++candidateIndex) {
      const Neighbor &neighbor = candidates[candidateIndex];
      if (workspaceEdgeBlocked(workspace, x, y, neighbor.x, neighbor.y))
        continue;

      const RoutePoint nextPoint{workspace.xCoordinates[neighbor.x],
                                 workspace.yCoordinates[neighbor.y]};
      double stepCost = std::abs(nextPoint.x - currentPoint.x) +
                        std::abs(nextPoint.y - currentPoint.y);
      if (direction != 0 && direction != neighbor.direction)
        stepCost += turnCost;
      stepCost +=
          workspaceWireCost(workspace, occupancy, x, y, neighbor.x, neighbor.y);

      const std::size_t nextState =
          encode(neighbor.x, neighbor.y, neighbor.direction);
      const double candidateDistance = distance[item.state] + stepCost;
      if (candidateDistance >= distance[nextState])
        continue;

      distance[nextState] = candidateDistance;
      predecessor[nextState] = item.state;
      queue.push(
          {candidateDistance + heuristic(neighbor.x, neighbor.y), nextState});
    }
  }

  if (goalState == std::numeric_limits<std::size_t>::max())
    return {};

  std::vector<RoutePoint> gridPoints;
  for (std::size_t state = goalState;
       state != std::numeric_limits<std::size_t>::max();
       state = predecessor[state]) {
    std::size_t x = 0, y = 0, direction = 0;
    decode(state, x, y, direction);
    (void)direction;
    gridPoints.push_back(
        {workspace.xCoordinates[x], workspace.yCoordinates[y]});
    if (state == startState)
      break;
  }

  if (gridPoints.empty() || !samePoint(gridPoints.back(), endpoints.escape))
    return {};

  std::reverse(gridPoints.begin(), gridPoints.end());
  RoutePath route;
  route.bus = connection.bus;
  route.component = connection.component;
  appendRoutePoint(route.points, endpoints.bus);
  for (const auto &point : gridPoints)
    appendRoutePoint(route.points, point);
  appendRoutePoint(route.points, endpoints.terminal);
  route.points = simplifyRoutePoints(route.points);
  route.score = distance[goalState];
  return route;
}

std::vector<std::vector<RouteBridge>>
findRouteBridges(const std::vector<RoutePath> &routes) {

  std::vector<std::vector<RouteBridge>> bridges(routes.size());

  for (std::size_t laterIndex = 0; laterIndex < routes.size(); ++laterIndex) {

    const auto laterSegments = routeSegments(routes[laterIndex].points);

    for (std::size_t earlierIndex = 0; earlierIndex < laterIndex;
         ++earlierIndex) {

      const auto earlierSegments = routeSegments(routes[earlierIndex].points);

      for (std::size_t laterSegmentIndex = 0;
           laterSegmentIndex < laterSegments.size(); ++laterSegmentIndex) {

        for (const auto &earlierSegment : earlierSegments) {

          RoutePoint crossing;

          if (!routeProperCrossing(laterSegments[laterSegmentIndex],
                                   earlierSegment, &crossing)) {

            continue;
          }

          const Bool horizontal =
              routeSegmentHorizontal(laterSegments[laterSegmentIndex]);

          bridges[laterIndex].push_back(
              {laterSegmentIndex, horizontal ? crossing.x : crossing.y});
        }
      }
    }
  }

  return bridges;
}

struct PendingRoute {
  std::size_t connectionIndex = 0;
  double port = 0.0;
};

void routeBusConnections(String &svg, std::vector<Connection> &connections) {

  // ---------------------------------------------------------------------------
  // This pass is deliberately independent of the selected topology layout.
  // LeftToRight, TopToBottom and Quadratic all receive exactly the same
  // collision checker and line router.
  // ---------------------------------------------------------------------------

  std::unordered_map<String, std::vector<std::size_t>> byBus;

  for (std::size_t index = 0; index < connections.size(); ++index) {

    byBus[connections[index].bus].push_back(index);
  }

  std::vector<PendingRoute> pendingRoutes;

  for (auto &busEntry : byBus) {

    auto &busIndices = busEntry.second;

    if (busIndices.empty()) {
      continue;
    }

    const BusbarOrientation busOrientation =
        connections[busIndices.front()].busOrientation;

    std::vector<std::size_t> negative;

    std::vector<std::size_t> positive;

    for (const std::size_t index : busIndices) {

      auto &connection = connections[index];

      connection.side = broadSideFor(
          connection.busBounds, connection.componentBounds, busOrientation);

      if (connection.side == BroadSide::Negative) {

        negative.push_back(index);

      } else {

        positive.push_back(index);
      }
    }

    auto prepareSide = [&](std::vector<std::size_t> &indices) {
      std::sort(indices.begin(), indices.end(),
                [&](const std::size_t first, const std::size_t second) {
                  return longCoordinate(connections[first].componentBounds,
                                        busOrientation) <
                         longCoordinate(connections[second].componentBounds,
                                        busOrientation);
                });

      const std::size_t count = indices.size();

      for (std::size_t position = 0; position < count; ++position) {

        const auto &connection = connections[indices[position]];

        const double fraction =
            static_cast<double>(position + 1) / static_cast<double>(count + 1);

        const double port =
            busLongMin(connection.busBounds, busOrientation) +
            fraction * busLongLength(connection.busBounds, busOrientation);

        pendingRoutes.push_back({indices[position], port});
      }
    };

    prepareSide(negative);

    prepareSide(positive);
  }

  const auto obstacles = routeObstacles(connections);

  std::vector<WorkspaceEndpointPair> workspaceEndpoints;
  workspaceEndpoints.reserve(pendingRoutes.size());
  for (const auto &pending : pendingRoutes) {
    const auto endpoints =
        routeEndpoints(connections[pending.connectionIndex], pending.port);
    workspaceEndpoints.push_back({endpoints.escape, endpoints.approach});
  }

  const RoutingWorkspace routingWorkspace =
      buildRoutingWorkspace(obstacles, workspaceEndpoints);

  RoutingOccupancy routingOccupancy = makeRoutingOccupancy(routingWorkspace);

  std::vector<RoutePath> routes;

  routes.reserve(pendingRoutes.size());

  std::vector<RouteSegment> occupiedSegments;

  occupiedSegments.reserve(pendingRoutes.size() * 4);

  UInt cheapRouteCount = 0;
  UInt gridRouteCount = 0;
  UInt fallbackRouteCount = 0;

  // Shorter connections are routed first. Long/difficult connections can then
  // choose outer lanes around the already occupied local geometry.
  std::sort(pendingRoutes.begin(), pendingRoutes.end(),
            [&](const PendingRoute &first, const PendingRoute &second) {
              const auto firstEndpoints = routeEndpoints(
                  connections[first.connectionIndex], first.port);

              const auto secondEndpoints = routeEndpoints(
                  connections[second.connectionIndex], second.port);

              const double firstDistance =
                  std::abs(firstEndpoints.terminal.x - firstEndpoints.bus.x) +
                  std::abs(firstEndpoints.terminal.y - firstEndpoints.bus.y);

              const double secondDistance =
                  std::abs(secondEndpoints.terminal.x - secondEndpoints.bus.x) +
                  std::abs(secondEndpoints.terminal.y - secondEndpoints.bus.y);

              return firstDistance < secondDistance;
            });

  for (const auto &pending : pendingRoutes) {

    const auto &connection = connections[pending.connectionIndex];

    // Fast path: most well-placed connections need only one or two Manhattan
    // bends. Avoid A* entirely when that route does not hit an object or wire.
    RoutePath route = cheapOrthogonalRoute(connection, pending.port, obstacles,
                                           occupiedSegments);

    if (!route.points.empty())
      ++cheapRouteCount;

    // Expensive path: invoke the cached visibility-grid/A* router only when the cheap
    // route collides. The A* search itself considers only a local corridor and
    // expands that corridor when necessary.
    if (route.points.empty()) {

      route = orthogonalGridRoute(connection, pending.port, routingWorkspace,
                                  routingOccupancy);
      if (!route.points.empty())
        ++gridRouteCount;
    }

    // Robust fallback for pathological geometry.
    if (route.points.empty()) {

      route = chooseBestRoute(connection, pending.port, obstacles, routes);
      ++fallbackRouteCount;
    }

    // HARD busbar rule: every connection MUST leave through a long side.
    if (!routeStartsAtLongBusSide(connection, route, pending.port)) {
      const RouteEndpoints endpoints = routeEndpoints(connection, pending.port);
      RoutePath safeRoute;
      safeRoute.bus = connection.bus;
      safeRoute.component = connection.component;
      appendRoutePoint(safeRoute.points, endpoints.bus);
      appendRoutePoint(safeRoute.points, endpoints.escape);
      if (connection.busOrientation == BusbarOrientation::Vertical) {
        // Keep the vertical run OUTSIDE the busbar; never use its top/bottom end.
        appendRoutePoint(safeRoute.points,
                         {endpoints.escape.x, endpoints.approach.y});
      } else {
        // Rotated equivalent: keep the horizontal run outside the busbar.
        appendRoutePoint(safeRoute.points,
                         {endpoints.approach.x, endpoints.escape.y});
      }
      appendRoutePoint(safeRoute.points, endpoints.approach);
      appendRoutePoint(safeRoute.points, endpoints.terminal);
      safeRoute.points = simplifyRoutePoints(safeRoute.points);
      route = std::move(safeRoute);
    }

    markRouteInOccupancy(routingWorkspace, routingOccupancy, route);

    const auto newSegments = routeSegments(route.points);

    occupiedSegments.insert(occupiedSegments.end(), newSegments.begin(),
                            newSegments.end());

    routes.push_back(std::move(route));
  }

  SPDLOG_INFO(
      "Topology router: connections={}, cheap={}, grid={}, fallback={}, "
      "workspace={}x{}.",
      pendingRoutes.size(), cheapRouteCount, gridRouteCount, fallbackRouteCount,
      routingWorkspace.xCount, routingWorkspace.yCount);

  // Any crossing that survived the route optimisation is treated as
  // unavoidable and receives a standard electrical drawing bridge/hump.
  const auto bridges = findRouteBridges(routes);

  for (std::size_t index = 0; index < routes.size(); ++index) {

    replaceEdgePath(svg, routes[index].bus, routes[index].component,
                    routePathWithBridges(routes[index], bridges[index]));
  }
}

// -----------------------------------------------------------------------------
// Ground symbols
// -----------------------------------------------------------------------------

GroundPlacement makeGroundPlacement(const std::filesystem::path &directory,
                                    const SvgBounds &component,
                                    const SvgBounds &bus,
                                    const SymbolOrientation orientation) {

  constexpr double groundSize = 10.0;
  constexpr double gap = 8.0;

  GroundPlacement placement;

  if (orientation == SymbolOrientation::Horizontal) {

    const Bool busLeft = centerX(bus) < centerX(component);

    if (busLeft) {
      placement.symbol = groundSymbol(directory, "left");

      placement.bounds.minX = component.maxX + gap;

      placement.bounds.maxX = placement.bounds.minX + groundSize;

      placement.bounds.minY = centerY(component) - 0.5 * groundSize;

      placement.bounds.maxY = placement.bounds.minY + groundSize;

      placement.lineX0 = component.maxX;
      placement.lineY0 = centerY(component);
      placement.lineX1 = placement.bounds.minX;
      placement.lineY1 = centerY(placement.bounds);

    } else {
      placement.symbol = groundSymbol(directory, "right");

      placement.bounds.maxX = component.minX - gap;

      placement.bounds.minX = placement.bounds.maxX - groundSize;

      placement.bounds.minY = centerY(component) - 0.5 * groundSize;

      placement.bounds.maxY = placement.bounds.minY + groundSize;

      placement.lineX0 = component.minX;
      placement.lineY0 = centerY(component);
      placement.lineX1 = placement.bounds.maxX;
      placement.lineY1 = centerY(placement.bounds);
    }

  } else {
    const Bool busAbove = centerY(bus) < centerY(component);

    if (busAbove) {
      placement.symbol = groundSymbol(directory, "top");

      placement.bounds.minY = component.maxY + gap;

      placement.bounds.maxY = placement.bounds.minY + groundSize;

      placement.bounds.minX = centerX(component) - 0.5 * groundSize;

      placement.bounds.maxX = placement.bounds.minX + groundSize;

      placement.lineX0 = centerX(component);
      placement.lineY0 = component.maxY;
      placement.lineX1 = centerX(placement.bounds);
      placement.lineY1 = placement.bounds.minY;

    } else {
      placement.symbol = groundSymbol(directory, "bottom");

      placement.bounds.maxY = component.minY - gap;

      placement.bounds.minY = placement.bounds.maxY - groundSize;

      placement.bounds.minX = centerX(component) - 0.5 * groundSize;

      placement.bounds.maxX = placement.bounds.minX + groundSize;

      placement.lineX0 = centerX(component);
      placement.lineY0 = component.minY;
      placement.lineX1 = centerX(placement.bounds);
      placement.lineY1 = placement.bounds.maxY;
    }
  }

  placement.bounds.valid = true;
  return placement;
}

void injectGroundConnection(String &svg, const GroundPlacement &placement) {

  const auto graphEnd = svg.rfind("</g>");

  if (graphEnd == String::npos)
    return;

  std::ostringstream path;

  path << "\n<path class=\"dpsim-ground-connection\""
       << " fill=\"none\""
       << " stroke=\"black\""
       << " stroke-width=\"0.25\""
       << " d=\"M" << formatReal(placement.lineX0) << ","
       << formatReal(placement.lineY0) << " L" << formatReal(placement.lineX1)
       << "," << formatReal(placement.lineY1) << "\"/>\n";

  svg.insert(graphEnd, path.str());
}

// -----------------------------------------------------------------------------
// Compact collision-aware label placement
// -----------------------------------------------------------------------------

enum class LabelSide { Below, Above, Right, Left };

struct Segment {
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
};

Bool overlaps(const SvgBounds &first, const SvgBounds &second,
              const double clearance = 0.0) {

  return first.valid && second.valid && first.maxX + clearance > second.minX &&
         first.minX - clearance < second.maxX &&
         first.maxY + clearance > second.minY &&
         first.minY - clearance < second.maxY;
}

Bool segmentHitsBox(const Segment &segment, SvgBounds box,
                    const double clearance) {

  box.minX -= clearance;
  box.maxX += clearance;
  box.minY -= clearance;
  box.maxY += clearance;

  // Orthogonal electrical paths dominate this renderer. A simple interval test
  // is therefore sufficient and intentionally conservative for label placement.
  const double minX = std::min(segment.x0, segment.x1);

  const double maxX = std::max(segment.x0, segment.x1);

  const double minY = std::min(segment.y0, segment.y1);

  const double maxY = std::max(segment.y0, segment.y1);

  if (std::abs(segment.y0 - segment.y1) < 1.0e-9) {

    return segment.y0 >= box.minY && segment.y0 <= box.maxY &&
           maxX >= box.minX && minX <= box.maxX;
  }

  if (std::abs(segment.x0 - segment.x1) < 1.0e-9) {

    return segment.x0 >= box.minX && segment.x0 <= box.maxX &&
           maxY >= box.minY && minY <= box.maxY;
  }

  // Conservative fallback for any diagonal/cubic approximation.
  return maxX >= box.minX && minX <= box.maxX && maxY >= box.minY &&
         minY <= box.maxY;
}

std::vector<Segment> electricalSegments(const String &svg) {

  std::vector<Segment> result;

  const std::regex pathRegex(R"regex(<path[^>]*\bd="([^"]+)"[^>]*/>)regex");

  const std::regex pointRegex(R"regex(([-+0-9.eE]+),([-+0-9.eE]+))regex");

  for (auto pathIterator =
           std::sregex_iterator(svg.begin(), svg.end(), pathRegex);
       pathIterator != std::sregex_iterator(); ++pathIterator) {

    const String data = (*pathIterator)[1].str();

    std::vector<std::pair<double, double>> points;

    for (auto pointIterator =
             std::sregex_iterator(data.begin(), data.end(), pointRegex);
         pointIterator != std::sregex_iterator(); ++pointIterator) {

      points.emplace_back(std::stod((*pointIterator)[1].str()),
                          std::stod((*pointIterator)[2].str()));
    }

    for (std::size_t index = 1; index < points.size(); ++index) {

      result.push_back({points[index - 1].first, points[index - 1].second,
                        points[index].first, points[index].second});
    }
  }

  return result;
}

double estimatedTextWidth(const String &text, const double fontSize) {

  const String plain =
      std::regex_replace(text, std::regex(R"regex(<[^>]*>)regex"), "");

  return std::max(fontSize,
                  static_cast<double>(plain.size()) * fontSize * 0.56);
}

Bool placeNodeLabel(String &svg, const String &nodeName,
                    const double primarySize, const double secondarySize,
                    const std::vector<LabelSide> &preference,
                    const std::vector<SvgBounds> &geometry,
                    const std::vector<Segment> &segments,
                    std::vector<SvgBounds> &placedLabels) {

  const String title = "<title>" + nodeName + "</title>";

  const auto titlePosition = svg.find(title);

  if (titlePosition == String::npos)
    return false;

  auto nodeEnd = svg.find("\n<!--", titlePosition);

  if (nodeEnd == String::npos)
    nodeEnd = svg.size();

  const String nodeBlock = svg.substr(titlePosition, nodeEnd - titlePosition);

  const SvgBounds nodeBounds = extractPolygonBounds(nodeBlock);

  if (!nodeBounds.valid)
    return false;

  const std::regex textRegex(R"regex(<text[^>]*>.*?</text>)regex");

  struct TextElement {
    std::size_t relativePosition = 0;
    std::size_t length = 0;
    String opening;
    String remainder;
    double size = 0.0;
  };

  std::vector<TextElement> elements;

  double labelWidth = 0.0;
  double labelHeight = 0.0;
  std::size_t index = 0;

  for (auto iterator =
           std::sregex_iterator(nodeBlock.begin(), nodeBlock.end(), textRegex);
       iterator != std::sregex_iterator(); ++iterator) {

    const String element = iterator->str();

    const auto openingEnd = element.find('>');

    const auto closingStart = element.rfind("</text>");

    if (openingEnd == String::npos || closingStart == String::npos)
      continue;

    const double size = index == 0 ? primarySize : secondarySize;

    const String visible =
        element.substr(openingEnd + 1, closingStart - openingEnd - 1);

    labelWidth = std::max(labelWidth, estimatedTextWidth(visible, size));

    labelHeight += size + (index == 0 ? 0.0 : 0.7);

    elements.push_back({static_cast<std::size_t>(iterator->position()),
                        static_cast<std::size_t>(iterator->length()),
                        element.substr(0, openingEnd + 1),
                        element.substr(openingEnd + 1), size});

    ++index;
  }

  if (elements.empty())
    return false;

  struct Candidate {
    SvgBounds bounds;
    LabelSide side = LabelSide::Below;
    double x = 0.0;
    double topY = 0.0;
    String anchor = "middle";
    double score = std::numeric_limits<double>::max();
  };

  Candidate best;

  for (const double gap :
       {2.2, 6.2, 10.2, 16.0, 24.0, 34.0, 46.0, 62.0, 82.0}) {
    for (const auto side : preference) {
      Candidate candidate;
      candidate.side = side;
      candidate.bounds.valid = true;

      const double x = centerX(nodeBounds);

      const double y = centerY(nodeBounds);

      switch (side) {
      case LabelSide::Below:
        candidate.bounds.minX = x - 0.5 * labelWidth;
        candidate.bounds.maxX = x + 0.5 * labelWidth;
        candidate.bounds.minY = nodeBounds.maxY + gap;
        candidate.bounds.maxY = candidate.bounds.minY + labelHeight;
        candidate.x = x;
        candidate.topY = candidate.bounds.minY;
        candidate.anchor = "middle";
        break;

      case LabelSide::Above:
        candidate.bounds.minX = x - 0.5 * labelWidth;
        candidate.bounds.maxX = x + 0.5 * labelWidth;
        candidate.bounds.maxY = nodeBounds.minY - gap;
        candidate.bounds.minY = candidate.bounds.maxY - labelHeight;
        candidate.x = x;
        candidate.topY = candidate.bounds.minY;
        candidate.anchor = "middle";
        break;

      case LabelSide::Right:
        candidate.bounds.minX = nodeBounds.maxX + gap;
        candidate.bounds.maxX = candidate.bounds.minX + labelWidth;
        candidate.bounds.minY = y - 0.5 * labelHeight;
        candidate.bounds.maxY = candidate.bounds.minY + labelHeight;
        candidate.x = candidate.bounds.minX;
        candidate.topY = candidate.bounds.minY;
        candidate.anchor = "start";
        break;

      case LabelSide::Left:
        candidate.bounds.maxX = nodeBounds.minX - gap;
        candidate.bounds.minX = candidate.bounds.maxX - labelWidth;
        candidate.bounds.minY = y - 0.5 * labelHeight;
        candidate.bounds.maxY = candidate.bounds.minY + labelHeight;
        candidate.x = candidate.bounds.maxX;
        candidate.topY = candidate.bounds.minY;
        candidate.anchor = "end";
        break;
      }

      double score = 0.0;

      for (const auto &box : geometry) {
        if (overlaps(candidate.bounds, box, 1.0))
          score += 1000.0;
      }

      for (const auto &box : placedLabels) {
        if (overlaps(candidate.bounds, box, 1.0))
          score += 1200.0;
      }

      for (const auto &segment : segments) {
        if (segmentHitsBox(segment, candidate.bounds, 0.8))
          score += 100.0;
      }

      candidate.score = score;

      if (score < best.score)
        best = candidate;

      if (score <= 1.0e-9)
        goto candidate_selected;
    }
  }

candidate_selected:

  struct Replacement {
    std::size_t position = 0;
    std::size_t length = 0;
    String value;
  };

  std::vector<Replacement> replacements;
  double yOffset = 0.0;

  for (const auto &element : elements) {
    String opening = element.opening;

    auto replaceOrInsert = [&opening](const std::regex &attribute,
                                      const String &value) {
      if (std::regex_search(opening, attribute)) {

        opening = std::regex_replace(opening, attribute, value);

      } else {
        opening.insert(opening.size() - 1, value);
      }
    };

    const double baseline = best.topY + yOffset + element.size;

    replaceOrInsert(std::regex(R"regex(\s+x="[^"]*")regex"),
                    " x=\"" + formatReal(best.x) + "\"");

    replaceOrInsert(std::regex(R"regex(\s+y="[^"]*")regex"),
                    " y=\"" + formatReal(baseline) + "\"");

    replaceOrInsert(std::regex(R"regex(\s+text-anchor="[^"]*")regex"),
                    " text-anchor=\"" + best.anchor + "\"");

    replaceOrInsert(std::regex(R"regex(\s+font-size="[^"]*")regex"),
                    " font-size=\"" + formatReal(element.size) + "\"");

    replaceOrInsert(std::regex(R"regex(\s+font-family="[^"]*")regex"),
                    " font-family=\"" + String(kFont) + "\"");

    replaceOrInsert(std::regex(R"regex(\s+font-weight="[^"]*")regex"),
                    " font-weight=\"normal\"");

    replacements.push_back({titlePosition + element.relativePosition,
                            element.length, opening + element.remainder});

    yOffset += element.size + 0.7;
  }

  for (auto iterator = replacements.rbegin(); iterator != replacements.rend();
       ++iterator) {

    svg.replace(iterator->position, iterator->length, iterator->value);
  }

  placedLabels.push_back(best.bounds);
  return true;
}

// -----------------------------------------------------------------------------
// Canvas crop
// -----------------------------------------------------------------------------

void includeBounds(SvgBounds &total, const SvgBounds &bounds) {

  if (!bounds.valid)
    return;

  if (!total.valid) {
    total = bounds;
    return;
  }

  total.minX = std::min(total.minX, bounds.minX);
  total.minY = std::min(total.minY, bounds.minY);
  total.maxX = std::max(total.maxX, bounds.maxX);
  total.maxY = std::max(total.maxY, bounds.maxY);
}

void fitCanvas(String &svg, const SvgBounds &content, const double margin) {

  if (!content.valid)
    return;

  const std::regex transformRegex(
      R"regex(<g[^>]*id="graph0"[^>]*transform="[^"]*translate\(([-+0-9.eE]+)[ ,]+([-+0-9.eE]+)\)[^"]*"[^>]*>)regex");

  std::smatch transformMatch;

  if (!std::regex_search(svg, transformMatch, transformRegex))
    return;

  const double translateX = std::stod(transformMatch[1].str());

  const double translateY = std::stod(transformMatch[2].str());

  // Keep a generous symmetric safety area around the final content.
  //
  // Text-width estimation is intentionally approximate and imported symbol SVGs
  // can contain strokes that slightly exceed the Graphviz placeholder bounds.
  // A larger fixed allowance is preferable to clipping publication figures.
  constexpr double canvasSafety = 28.0;

  const double minX = content.minX + translateX - margin - canvasSafety;

  const double minY = content.minY + translateY - margin - canvasSafety;

  const double maxX = content.maxX + translateX + margin + canvasSafety;

  const double maxY = content.maxY + translateY + margin + canvasSafety;

  const double width = std::max(1.0, maxX - minX);

  const double height = std::max(1.0, maxY - minY);

  const std::regex svgRegex(R"regex(<svg([^>]*)>)regex");

  std::smatch svgMatch;

  if (!std::regex_search(svg, svgMatch, svgRegex))
    return;

  String opening = svgMatch[0].str();

  const String viewBox = "viewBox=\"" + formatReal(minX) + " " +
                         formatReal(minY) + " " + formatReal(width) + " " +
                         formatReal(height) + "\"";

  opening = std::regex_replace(
      opening,
      std::regex(
          R"regex(viewBox="[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+")regex"),
      viewBox);

  opening = std::regex_replace(opening,
                               std::regex(R"regex(width="[-+0-9.eE]+pt")regex"),
                               "width=\"" + formatReal(width) + "pt\"");

  opening = std::regex_replace(
      opening, std::regex(R"regex(height="[-+0-9.eE]+pt")regex"),
      "height=\"" + formatReal(height) + "pt\"");

  svg.replace(static_cast<std::size_t>(svgMatch.position()),
              static_cast<std::size_t>(svgMatch.length()), opening);

  // Graphviz's original white background polygon still has the ORIGINAL graph
  // dimensions. Once topology post-processing moves nodes outside those bounds,
  // parts of the enlarged viewBox can otherwise lie outside the white area.
  //
  // Insert a root-level white rectangle covering the complete FINAL viewBox.
  // It is inserted before graph0, so it can never cover electrical content.
  const auto rootSvgStart = svg.find("<svg");

  const auto rootSvgEnd =
      rootSvgStart == String::npos ? String::npos : svg.find('>', rootSvgStart);

  if (rootSvgEnd != String::npos) {

    const String whiteBackground =
        "\n<rect class=\"dpsim-canvas-background\""
        " x=\"" +
        formatReal(minX) + "\" y=\"" + formatReal(minY) + "\" width=\"" +
        formatReal(width) + "\" height=\"" + formatReal(height) +
        "\" fill=\"white\" stroke=\"none\"/>\n";

    svg.insert(rootSvgEnd + 1, whiteBackground);
  }
}

} // namespace

double centerX(const SvgBounds &bounds) {
  return 0.5 * (bounds.minX + bounds.maxX);
}

double centerY(const SvgBounds &bounds) {
  return 0.5 * (bounds.minY + bounds.maxY);
}

Bool extractNodeBounds(const String &svg, const String &graphNodeName,
                       SvgBounds &bounds) {

  const String title = "<title>" + graphNodeName + "</title>";

  const auto titlePosition = svg.find(title);

  if (titlePosition == String::npos)
    return false;

  auto nodeEnd = svg.find("\n<!--", titlePosition);

  if (nodeEnd == String::npos)
    nodeEnd = svg.size();

  bounds =
      extractPolygonBounds(svg.substr(titlePosition, nodeEnd - titlePosition));

  return bounds.valid;
}

Bool replaceNodeBounds(String &svg, const String &graphNodeName,
                       const SvgBounds &bounds) {

  const String title = "<title>" + graphNodeName + "</title>";

  const auto titlePosition = svg.find(title);

  if (titlePosition == String::npos)
    return false;

  auto nodeEnd = svg.find("\n<!--", titlePosition);

  if (nodeEnd == String::npos)
    nodeEnd = svg.size();

  String block = svg.substr(titlePosition, nodeEnd - titlePosition);

  const std::regex polygonRegex(R"regex(points="([^"]+)")regex");

  std::smatch match;

  if (!std::regex_search(block, match, polygonRegex))
    return false;

  const String replacement =
      "points=\"" + formatReal(bounds.minX) + "," + formatReal(bounds.minY) +
      " " + formatReal(bounds.maxX) + "," + formatReal(bounds.minY) + " " +
      formatReal(bounds.maxX) + "," + formatReal(bounds.maxY) + " " +
      formatReal(bounds.minX) + "," + formatReal(bounds.maxY) + " " +
      formatReal(bounds.minX) + "," + formatReal(bounds.minY) + "\"";

  block.replace(static_cast<std::size_t>(match.position()),
                static_cast<std::size_t>(match.length()), replacement);

  svg.replace(titlePosition, nodeEnd - titlePosition, block);

  return true;
}

Bool moveNodeCenter(String &svg, const String &graphNodeName, const double x,
                    const double y) {

  SvgBounds bounds;

  if (!extractNodeBounds(svg, graphNodeName, bounds))
    return false;

  const double dx = x - centerX(bounds);

  const double dy = y - centerY(bounds);

  bounds.minX += dx;
  bounds.maxX += dx;
  bounds.minY += dy;
  bounds.maxY += dy;

  return replaceNodeBounds(svg, graphNodeName, bounds);
}

Bool isLineLike(TopologicalPowerComp &component) {

  const String type = component.type();

  return contains(type, "PiLine") || contains(type, "RxLine") ||
         contains(type, "RXLine") || contains(type, "TransmissionLine") ||
         contains(type, "Cable");
}

String postProcessSvg(String svg, const SystemTopology &system,
                      const std::filesystem::path &symbolDirectory,
                      const SystemTopologyRenderer::Options &options,
                      const LayoutResult &layoutResult) {

  const auto postProcessStart = std::chrono::steady_clock::now();

  std::unordered_map<String, SymbolOrientation> componentOrientations;

  std::unordered_map<String, Bool> componentHasGround;

  std::unordered_map<String, Bool> componentSingleBus;

  std::unordered_map<String, String> componentPhysicalBus;

  std::unordered_map<String, UInt> groundedLeafCount;

  // ---------------------------------------------------------------------------
  // Component orientation + local leaf placement
  // ---------------------------------------------------------------------------

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    if (RenderStyle::collapseLines && isLineLike(*component))
      continue;

    const String componentNodeName = "component_" + component->uid();

    SvgBounds componentBounds;

    if (!extractNodeBounds(svg, componentNodeName, componentBounds))
      continue;

    std::vector<SvgBounds> busBounds;
    Bool hasGround = false;
    String physicalBus;

    for (const auto &terminal : component->topologicalTerminals()) {
      if (!terminal)
        continue;

      const auto node = terminal->topologicalNodes();

      if (!node)
        continue;

      if (node->isGround()) {
        hasGround = true;
        continue;
      }

      const String busName = "bus_" + node->uid();

      SvgBounds bounds;

      if (extractNodeBounds(svg, busName, bounds)) {

        busBounds.push_back(bounds);

        if (physicalBus.empty())
          physicalBus = busName;
      }
    }

    const Bool groundedLeaf = hasGround && busBounds.size() == 1;

    const Bool singleBus = busBounds.size() == 1;

    SymbolOrientation orientation = SymbolOrientation::Horizontal;

    if (groundedLeaf && !physicalBus.empty()) {

      auto directionIterator = layoutResult.leafDirections.find(physicalBus);

      Direction direction;

      if (directionIterator != layoutResult.leafDirections.end()) {

        direction = directionIterator->second;

      } else {
        direction =
            options.layout == SystemTopologyRenderer::Layout::LeftToRight
                ? Direction{0.0, -1.0}
                : Direction{1.0, 0.0};
      }

      orientation = std::abs(direction.x) >= std::abs(direction.y)
                        ? SymbolOrientation::Horizontal
                        : SymbolOrientation::Vertical;

    } else if (busBounds.size() >= 2) {
      double largestDistance = -1.0;
      double selectedDx = 0.0;
      double selectedDy = 0.0;

      for (std::size_t first = 0; first < busBounds.size(); ++first) {

        for (std::size_t second = first + 1; second < busBounds.size();
             ++second) {

          const double dx =
              centerX(busBounds[second]) - centerX(busBounds[first]);

          const double dy =
              centerY(busBounds[second]) - centerY(busBounds[first]);

          const double distance = dx * dx + dy * dy;

          if (distance > largestDistance) {
            largestDistance = distance;
            selectedDx = dx;
            selectedDy = dy;
          }
        }
      }

      if (std::abs(selectedDy) > std::abs(selectedDx))
        orientation = SymbolOrientation::Vertical;

    } else if (busBounds.size() == 1) {
      const double dx = centerX(busBounds.front()) - centerX(componentBounds);

      const double dy = centerY(busBounds.front()) - centerY(componentBounds);

      if (std::abs(dy) > std::abs(dx))
        orientation = SymbolOrientation::Vertical;
    }

    componentOrientations[componentNodeName] = orientation;

    componentHasGround[componentNodeName] = hasGround;

    componentSingleBus[componentNodeName] = singleBus;

    if (!physicalBus.empty())
      componentPhysicalBus[componentNodeName] = physicalBus;

    SvgBounds resized = resizedBounds(
        componentBounds, orientedDimensions(component->type(), orientation));

    if (groundedLeaf && !physicalBus.empty()) {

      SvgBounds bus;

      if (extractNodeBounds(svg, physicalBus, bus)) {

        Direction direction;

        const auto iterator = layoutResult.leafDirections.find(physicalBus);

        if (iterator != layoutResult.leafDirections.end()) {

          direction = cardinal(iterator->second);

        } else {
          direction =
              options.layout == SystemTopologyRenderer::Layout::LeftToRight
                  ? Direction{0.0, -1.0}
                  : Direction{1.0, 0.0};
        }

        UInt &index = groundedLeafCount[physicalBus];

        // First leaf follows preferred free direction; additional leaves
        // alternate to the opposite side, then use parallel lanes.
        const double sign = index % 2 == 0 ? 1.0 : -1.0;

        Direction actual{direction.x * sign, direction.y * sign};

        const UInt layer = index / 2;

        const double extra = static_cast<double>(layer) * 34.0;

        constexpr double gap = 32.0;

        double targetX = centerX(bus);

        double targetY = centerY(bus);

        if (std::abs(actual.x) >= std::abs(actual.y)) {

          targetX +=
              actual.x * (0.5 * bus.width() + gap + 0.5 * resized.width());

          targetY += extra;

        } else {
          targetY +=
              actual.y * (0.5 * bus.height() + gap + 0.5 * resized.height());

          targetX += extra;
        }

        const double dx = targetX - centerX(resized);

        const double dy = targetY - centerY(resized);

        resized.minX += dx;
        resized.maxX += dx;
        resized.minY += dy;
        resized.maxY += dy;

        ++index;
      }
    }

    replaceNodeBounds(svg, componentNodeName, resized);
  }

  // ---------------------------------------------------------------------------
  // Busbar orientation
  // ---------------------------------------------------------------------------

  std::unordered_map<String, BusbarOrientation> busOrientations;

  for (const auto &node : system.mNodes) {
    if (!node || node->isGround())
      continue;

    const String busName = "bus_" + node->uid();

    BusbarOrientation orientation =
        options.layout == SystemTopologyRenderer::Layout::LeftToRight
            ? BusbarOrientation::Vertical
            : BusbarOrientation::Horizontal;

    if (layoutResult.circular) {

      // For a closed ring, the busbar orientation is derived from the OUTWARD
      // direction of the cycle.
      //
      // This is the robust rule for corner buses:
      //
      //   outward left/right -> VERTICAL busbar
      //   outward up/down    -> HORIZONTAL busbar
      //
      // The long side therefore faces the cycle interior/exterior and every
      // external branch/load leaves from a broad side. This restores the
      // orientation that worked correctly for IEEE-9 BUS4/BUS5/BUS6/BUS7/
      // BUS8/BUS9 before the last refactor.
      const auto outwardIterator = layoutResult.leafDirections.find(busName);

      if (outwardIterator != layoutResult.leafDirections.end()) {

        const Direction outward = outwardIterator->second;

        orientation = std::abs(outward.x) >= std::abs(outward.y)
                          ? BusbarOrientation::Vertical
                          : BusbarOrientation::Horizontal;
      }

    } else {

      // Non-circular grid:
      // orient the busbar perpendicular to the longest-path/backbone direction.
      const auto directionIterator = layoutResult.mainDirections.find(busName);

      if (directionIterator != layoutResult.mainDirections.end()) {

        const Direction direction = directionIterator->second;

        // Main electrical path horizontal -> vertical busbar.
        // Main electrical path vertical   -> horizontal busbar.
        orientation = std::abs(direction.x) >= std::abs(direction.y)
                          ? BusbarOrientation::Vertical
                          : BusbarOrientation::Horizontal;
      }
    }

    busOrientations[busName] = orientation;

    SvgBounds bounds;

    if (!extractNodeBounds(svg, busName, bounds))
      continue;

    const SymbolDimensions dimensions =
        orientation == BusbarOrientation::Vertical
            ? SymbolDimensions{0.075, 0.52}
            : SymbolDimensions{0.52, 0.075};

    replaceNodeBounds(svg, busName, resizedBounds(bounds, dimensions));
  }

  // ---------------------------------------------------------------------------
  // Final object-collision validation pass
  //
  // Applies to LeftToRight, TopToBottom and Quadratic alike. The topology-aware
  // 2-core/block layout has already decided the graph geometry. This pass only
  // creates additional space and separates local component collisions; it does
  // not independently rearrange overlapping backbone buses.
  // ---------------------------------------------------------------------------

  const auto overlapStart = std::chrono::steady_clock::now();

  const UInt remainingObjectOverlaps = resolveObjectOverlaps(svg, system);

  const auto overlapEnd = std::chrono::steady_clock::now();

  if (remainingObjectOverlaps > 0) {

    SPDLOG_WARN("Topology renderer: {} object overlaps remain after collision "
                "resolution.",
                remainingObjectOverlaps);
  }

  // ---------------------------------------------------------------------------
  // Build and route every bus-component connection.
  // The router only uses the LONG busbar sides.
  // ---------------------------------------------------------------------------

  std::vector<Connection> connections;

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    if (RenderStyle::collapseLines && isLineLike(*component))
      continue;

    const String componentName = "component_" + component->uid();

    SvgBounds componentBounds;

    if (!extractNodeBounds(svg, componentName, componentBounds))
      continue;

    for (const auto &terminal : component->topologicalTerminals()) {
      if (!terminal)
        continue;

      const auto node = terminal->topologicalNodes();

      if (!node || node->isGround())
        continue;

      const String busName = "bus_" + node->uid();

      SvgBounds busBounds;

      if (!extractNodeBounds(svg, busName, busBounds))
        continue;

      const auto busOrientationIterator = busOrientations.find(busName);

      const auto componentOrientationIterator =
          componentOrientations.find(componentName);

      if (busOrientationIterator == busOrientations.end() ||
          componentOrientationIterator == componentOrientations.end()) {

        SPDLOG_WARN("Topology renderer: skipping connection {} <-> {} because "
                    "orientation geometry is unavailable.",
                    busName, componentName);

        continue;
      }

      Connection connection;
      connection.bus = busName;
      connection.component = componentName;
      connection.busBounds = busBounds;
      connection.componentBounds = componentBounds;
      connection.busOrientation = busOrientationIterator->second;
      connection.componentOrientation = componentOrientationIterator->second;
      connection.singleBusComponent = componentSingleBus[componentName];
      connection.groundedLeaf = componentHasGround[componentName];

      connections.push_back(std::move(connection));
    }
  }

  const auto routingStart = std::chrono::steady_clock::now();

  routeBusConnections(svg, connections);

  const auto routingEnd = std::chrono::steady_clock::now();

  // ---------------------------------------------------------------------------
  // Ground placement after bus routing / possible leaf alignment
  // ---------------------------------------------------------------------------

  std::vector<GroundPlacement> grounds;

  for (const auto &entry : componentHasGround) {
    if (!entry.second)
      continue;

    const String &componentName = entry.first;

    const auto busIterator = componentPhysicalBus.find(componentName);

    if (busIterator == componentPhysicalBus.end())
      continue;

    SvgBounds componentBounds;
    SvgBounds busBounds;

    if (!extractNodeBounds(svg, componentName, componentBounds) ||
        !extractNodeBounds(svg, busIterator->second, busBounds))
      continue;

    const auto orientationIterator = componentOrientations.find(componentName);

    if (orientationIterator == componentOrientations.end()) {

      SPDLOG_WARN(
          "Topology renderer: cannot place ground for {} because component "
          "orientation is unavailable.",
          componentName);

      continue;
    }

    grounds.push_back(makeGroundPlacement(symbolDirectory, componentBounds,
                                          busBounds,
                                          orientationIterator->second));
  }

  // Keep ground symbols clear of all buses/components as well.
  {
    std::vector<SvgBounds> groundObstacles;

    for (const auto &node : system.mNodes) {

      if (!node || node->isGround()) {

        continue;
      }

      SvgBounds bounds;

      if (extractNodeBounds(svg, "bus_" + node->uid(), bounds)) {

        groundObstacles.push_back(bounds);
      }
    }

    for (const auto &identifiedObject : system.mComponents) {

      if (!identifiedObject) {
        continue;
      }

      auto component =
          std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

      if (!component) {
        continue;
      }

      if (RenderStyle::collapseLines && isLineLike(*component)) {

        continue;
      }

      SvgBounds bounds;

      if (extractNodeBounds(svg, "component_" + component->uid(), bounds)) {

        groundObstacles.push_back(bounds);
      }
    }

    resolveGroundOverlaps(grounds, groundObstacles);
  }

  // ---------------------------------------------------------------------------
  // Label collision model
  // ---------------------------------------------------------------------------

  std::vector<SvgBounds> geometry;

  for (const auto &node : system.mNodes) {
    if (!node || node->isGround())
      continue;

    SvgBounds bounds;

    if (extractNodeBounds(svg, "bus_" + node->uid(), bounds))
      geometry.push_back(bounds);
  }

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    if (RenderStyle::collapseLines && isLineLike(*component))
      continue;

    SvgBounds bounds;

    if (extractNodeBounds(svg, "component_" + component->uid(), bounds))
      geometry.push_back(bounds);
  }

  for (const auto &ground : grounds)
    geometry.push_back(ground.bounds);

  std::vector<Segment> segments = electricalSegments(svg);

  for (const auto &ground : grounds) {
    segments.push_back(
        {ground.lineX0, ground.lineY0, ground.lineX1, ground.lineY1});
  }

  std::vector<SvgBounds> labelBounds;

  std::vector<LabelSide> componentPreference;
  std::vector<LabelSide> busPreference;

  if (options.layout == SystemTopologyRenderer::Layout::LeftToRight) {

    componentPreference = {LabelSide::Below, LabelSide::Above, LabelSide::Right,
                           LabelSide::Left};

    busPreference = componentPreference;

  } else {
    componentPreference = {LabelSide::Right, LabelSide::Left, LabelSide::Below,
                           LabelSide::Above};

    busPreference = {LabelSide::Below, LabelSide::Right, LabelSide::Left,
                     LabelSide::Above};
  }

  const auto labelStart = std::chrono::steady_clock::now();

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    if (RenderStyle::collapseLines && isLineLike(*component))
      continue;

    placeNodeLabel(svg, "component_" + component->uid(), kComponentLabelSize,
                   kTypeLabelSize, componentPreference, geometry, segments,
                   labelBounds);
  }

  for (const auto &node : system.mNodes) {
    if (!node || node->isGround())
      continue;

    placeNodeLabel(svg, "bus_" + node->uid(), kBusLabelSize, kBusLabelSize,
                   busPreference, geometry, segments, labelBounds);
  }

  const auto labelEnd = std::chrono::steady_clock::now();

  // ---------------------------------------------------------------------------
  // Content bounds before symbol injection
  // ---------------------------------------------------------------------------

  SvgBounds content;

  for (const auto &box : geometry)
    includeBounds(content, box);

  for (const auto &box : labelBounds)
    includeBounds(content, box);

  // Include every routed electrical segment in the final crop as well.
  // This prevents one-corner connections from being clipped after buses or
  // components have been moved by the topology-aware post-processing.
  for (const auto &segment : segments) {

    SvgBounds segmentBounds;

    segmentBounds.minX = std::min(segment.x0, segment.x1);

    segmentBounds.maxX = std::max(segment.x0, segment.x1);

    segmentBounds.minY = std::min(segment.y0, segment.y1);

    segmentBounds.maxY = std::max(segment.y0, segment.y1);

    // Zero-width / zero-height line bounds are still real content. Give them a
    // tiny thickness so includeBounds() treats them as valid.
    constexpr double linePadding = 0.5;

    if (segmentBounds.maxX - segmentBounds.minX < 1.0e-9) {

      segmentBounds.minX -= linePadding;

      segmentBounds.maxX += linePadding;
    }

    if (segmentBounds.maxY - segmentBounds.minY < 1.0e-9) {

      segmentBounds.minY -= linePadding;

      segmentBounds.maxY += linePadding;
    }

    segmentBounds.valid = true;

    includeBounds(content, segmentBounds);
  }

  // ---------------------------------------------------------------------------
  // Busbars
  // ---------------------------------------------------------------------------

  UInt embeddedSymbolCount = 0;

  if (RenderStyle::useSymbols) {
    for (const auto &node : system.mNodes) {
      if (!node || node->isGround())
        continue;

      const String busName = "bus_" + node->uid();

      const auto orientationIterator = busOrientations.find(busName);

      if (orientationIterator == busOrientations.end()) {

        SPDLOG_WARN(
            "Topology renderer: no orientation geometry for bus {}; keeping "
            "the Graphviz fallback representation.",
            busName);

        continue;
      }

      const auto symbol =
          busbarSymbol(symbolDirectory, orientationIterator->second);

      if (symbol.empty())
        continue;

      if (injectSymbolIntoNode(svg, busName, symbol, true))
        ++embeddedSymbolCount;
    }
  }

  // ---------------------------------------------------------------------------
  // Components
  // ---------------------------------------------------------------------------

  if (RenderStyle::useSymbols) {
    for (const auto &identifiedObject : system.mComponents) {
      if (!identifiedObject)
        continue;

      auto component =
          std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

      if (!component)
        continue;

      if (RenderStyle::collapseLines && isLineLike(*component))
        continue;

      const String componentName = "component_" + component->uid();

      const auto orientationIterator =
          componentOrientations.find(componentName);

      if (orientationIterator == componentOrientations.end()) {

        // Leave Graphviz's generic box in place. This can only happen when
        // Graphviz did not expose editable geometry for the component.
        SPDLOG_WARN(
            "Topology renderer: no orientation geometry for {}; keeping the "
            "Graphviz fallback representation.",
            componentName);

        continue;
      }

      const auto symbol = componentSymbol(symbolDirectory, *component,
                                          orientationIterator->second);

      if (symbol.empty())
        continue;

      if (injectSymbolIntoNode(svg, componentName, symbol))
        ++embeddedSymbolCount;
    }
  }

  // ---------------------------------------------------------------------------
  // Grounds
  // ---------------------------------------------------------------------------

  if (RenderStyle::useSymbols) {
    for (const auto &ground : grounds) {
      injectGroundConnection(svg, ground);

      if (injectFreeSymbol(svg, ground.symbol, ground.bounds, "dpsim-ground"))
        ++embeddedSymbolCount;
    }
  }

  fitCanvas(svg, content, 18.0);

  const auto postProcessEnd = std::chrono::steady_clock::now();

  const auto overlapMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(overlapEnd -
                                                            overlapStart)
          .count();

  const auto routingMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(routingEnd -
                                                            routingStart)
          .count();

  const auto labelMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(labelEnd -
                                                            labelStart)
          .count();

  const auto postProcessMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(postProcessEnd -
                                                            postProcessStart)
          .count();

  SPDLOG_INFO("Topology renderer: symbols={}, overlap={} ms, routing={} ms, "
              "labels={} ms, post-process={} ms.",
              embeddedSymbolCount, overlapMilliseconds, routingMilliseconds,
              labelMilliseconds, postProcessMilliseconds);

  return svg;
}

} // namespace CPS::TopologyRenderDetail

// =============================================================================
// SystemTopologyRenderer public implementation and Graphviz graph creation
// =============================================================================

namespace {

constexpr const char *kFont = "Latin Modern Roman";

constexpr double kComponentLabelSize = 5.5;
constexpr double kComponentTypeSize = 4.25;
constexpr double kBusLabelSize = 5.0;
constexpr double kLineLabelSize = 5.0;

Bool contains(const String &text, const String &needle) {

  return text.find(needle) != String::npos;
}

String formatReal(const double value) {

  std::ostringstream stream;
  stream << value;
  return stream.str();
}

String escapeLabel(const String &text) {

  String result;
  result.reserve(text.size());

  for (const char character : text) {
    switch (character) {
    case '&':
      result += "&amp;";
      break;
    case '<':
      result += "&lt;";
      break;
    case '>':
      result += "&gt;";
      break;
    case '"':
      result += "&quot;";
      break;
    default:
      result += character;
      break;
    }
  }

  return result;
}

String tooltipLowerCopy(String value) {

  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char character) {
                   return static_cast<char>(std::tolower(character));
                 });

  return value;
}

String trimTooltipText(String value) {

  const auto first = value.find_first_not_of(" \t\r\n");

  if (first == String::npos) {

    return {};
  }

  const auto last = value.find_last_not_of(" \t\r\n");

  return value.substr(first, last - first + 1);
}

Bool tooltipNameContains(const String &lowerName, const String &needle) {

  return lowerName.find(needle) != String::npos;
}

Bool isSetpointOrRatingAttribute(const String &attributeName) {

  const String name = tooltipLowerCopy(attributeName);

  // Configuration values that are especially useful in topology inspection.
  // Check these BEFORE runtime filtering so names such as `initial_voltage`,
  // `P_ref`, `Q_set`, `V_target`, etc. are never mistaken for measurements.
  const std::vector<String> configurationTokens{
      "_ref",   "ref_",  "reference", "_set",    "set_",    "setpoint",
      "target", "rated", "nominal",   "base",    "initial", "init_",
      "_init",  "limit", "maximum",   "minimum", "_max",    "_min"};

  for (const auto &token : configurationTokens) {

    if (tooltipNameContains(name, token)) {

      return true;
    }
  }

  // Common compact DPsim naming conventions.
  static const std::unordered_set<String> exactConfigurationNames{
      "p_ref",  "q_ref",  "s_ref",  "v_ref",  "i_ref",  "f_ref",
      "p_set",  "q_set",  "s_set",  "v_set",  "i_set",  "f_set",
      "pref",   "qref",   "sref",   "vref",   "iref",   "fref",
      "pset",   "qset",   "sset",   "vset",   "iset",   "fset",
      "p0",     "q0",     "s0",     "v0",     "i0",     "f0",
      "p_init", "q_init", "s_init", "v_init", "i_init", "f_init"};

  return exactConfigurationNames.find(name) != exactConfigurationNames.end();
}

Bool isRuntimeTooltipAttribute(const String &attributeName) {

  const String name = tooltipLowerCopy(attributeName);

  if (isSetpointOrRatingAttribute(attributeName)) {

    return false;
  }

  // Generic metadata / solver bookkeeping.
  static const std::unordered_set<String> exactRuntimeNames{
      "name",     "uid",    "left_vector", "right_vector", "state",  "states",
      "solution", "rhs",    "time",        "step",         "i",      "v",
      "p",        "q",      "s",           "i_intf",       "v_intf", "p_intf",
      "q_intf",   "p_elec", "q_elec",      "p_mech",       "t_e",    "delta_r",
      "w_r",      "theta",  "omega",       "frequency_pu"};

  if (exactRuntimeNames.find(name) != exactRuntimeNames.end()) {

    return true;
  }

  // Interface quantities and MNA/internal solver state are runtime results,
  // not configured component parameters.
  if (tooltipNameContains(name, "_intf") ||
      tooltipNameContains(name, "interface_current") ||
      tooltipNameContains(name, "interface_voltage") ||
      tooltipNameContains(name, "right_vector") ||
      tooltipNameContains(name, "left_vector") ||
      tooltipNameContains(name, "mna_") ||
      tooltipNameContains(name, "history") ||
      tooltipNameContains(name, "residual") ||
      tooltipNameContains(name, "derivative")) {

    return true;
  }

  // Common dynamic/result names. Keep explicit references/setpoints because
  // those are configured parameters and are handled separately below.
  if ((tooltipNameContains(name, "current") ||
       tooltipNameContains(name, "voltage") ||
       tooltipNameContains(name, "power") ||
       tooltipNameContains(name, "frequency") ||
       tooltipNameContains(name, "angle")) &&
      !tooltipNameContains(name, "ref") && !tooltipNameContains(name, "set") &&
      !tooltipNameContains(name, "rated") &&
      !tooltipNameContains(name, "nominal") &&
      !tooltipNameContains(name, "base") && !tooltipNameContains(name, "max") &&
      !tooltipNameContains(name, "min") &&
      !tooltipNameContains(name, "limit")) {

    return true;
  }

  return false;
}

Bool isLikelyConfiguredParameter(const String &attributeName) {

  const String name = tooltipLowerCopy(attributeName);

  if (isRuntimeTooltipAttribute(attributeName)) {

    return false;
  }

  if (isSetpointOrRatingAttribute(attributeName)) {

    return true;
  }

  // Explicit parameter semantics.
  const std::vector<String> parameterTokens{
      "resistance",  "inductance", "capacitance",   "conductance",  "reactance",
      "susceptance", "impedance",  "admittance",    "rated",        "nominal",
      "base",        "ratio",      "ref",           "set",          "limit",
      "maximum",     "minimum",    "max",           "min",          "gain",
      "droop",       "inertia",    "time_constant", "timeconstant", "tau",
      "frequency",   "voltage",    "current",       "power",        "phase",
      "pole",        "snubber"};

  for (const auto &token : parameterTokens) {

    if (tooltipNameContains(name, token)) {

      return true;
    }
  }

  // DPsim parameter attributes commonly use compact electrical/control names:
  // R_*, L_*, C_*, G_*, X_*, B_*, Kp/Ki, T*, H, etc.
  if (name.size() >= 2) {

    const char first = name.front();

    if ((first == 'r' || first == 'l' || first == 'c' || first == 'g' ||
         first == 'x' || first == 'b' || first == 'k' || first == 't') &&
        (name[1] == '_' || std::isdigit(static_cast<unsigned char>(name[1])) ||
         std::isalpha(static_cast<unsigned char>(name[1])))) {

      return true;
    }
  }

  if (name == "h" || name == "rs" || name == "ld" || name == "lq" ||
      name == "ll" || name == "xd" || name == "xq" || name == "kp" ||
      name == "ki" || name == "kd") {

    return true;
  }

  return false;
}

std::vector<double> tooltipNumericValues(const String &text) {

  // AttributeBase::toString() often prints Eigen vectors/matrices as rows of
  // whitespace-separated values. Extract their scalar entries so balanced
  // three-phase parameter matrices can be shown as one representative value.
  static const std::regex numberRegex(
      R"regex([-+]?(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[eE][-+]?[0-9]+)?)regex");

  std::vector<double> values;

  for (auto iterator =
           std::sregex_iterator(text.begin(), text.end(), numberRegex);
       iterator != std::sregex_iterator(); ++iterator) {

    try {
      values.push_back(std::stod(iterator->str()));
    } catch (...) {
      // Keep tooltip generation non-fatal. If parsing fails, the raw attribute
      // text is used below.
    }
  }

  return values;
}

String formatTooltipScalar(const double value) {

  std::ostringstream stream;

  stream << std::setprecision(7) << value;

  return stream.str();
}

Bool isCompactElectricalMatrixParameter(const String &attributeName) {

  const String name = tooltipLowerCopy(attributeName);

  return name == "r" || name == "l" || name == "c" || name == "g" ||
         name == "x" || name == "b" || name.rfind("r_", 0) == 0 ||
         name.rfind("l_", 0) == 0 || name.rfind("c_", 0) == 0 ||
         name.rfind("g_", 0) == 0 || name.rfind("x_", 0) == 0 ||
         name.rfind("b_", 0) == 0 || tooltipNameContains(name, "resistance") ||
         tooltipNameContains(name, "inductance") ||
         tooltipNameContains(name, "capacitance") ||
         tooltipNameContains(name, "conductance") ||
         tooltipNameContains(name, "reactance") ||
         tooltipNameContains(name, "susceptance");
}

String compactTooltipValue(const String &attributeName,
                           const String &rawValue) {

  const String trimmed = trimTooltipText(rawValue);

  if (trimmed.empty()) {
    return trimmed;
  }

  const Bool looksLikeMatrix =
      trimmed.find('\n') != String::npos || trimmed.find('\r') != String::npos;

  if (!looksLikeMatrix && !isCompactElectricalMatrixParameter(attributeName)) {

    return trimmed;
  }

  const auto values = tooltipNumericValues(trimmed);

  if (values.empty()) {
    return trimmed;
  }

  // For balanced three-phase R/L/C/G/X/B matrices, one representative
  // non-zero element is exactly what is useful in a topology tooltip.
  if (isCompactElectricalMatrixParameter(attributeName)) {

    constexpr double zeroTolerance = 1.0e-18;

    for (const double value : values) {

      if (std::abs(value) > zeroTolerance) {

        return formatTooltipScalar(value);
      }
    }

    return "0";
  }

  // For another multi-line numeric attribute, compact only if every non-zero
  // entry is effectively the same. Otherwise preserve the original value,
  // because reducing a genuinely non-uniform parameter matrix would hide data.
  Bool haveRepresentative = false;

  double representative = 0.0;

  for (const double value : values) {

    if (std::abs(value) < 1.0e-18) {

      continue;
    }

    if (!haveRepresentative) {

      representative = value;

      haveRepresentative = true;

      continue;
    }

    const double scale =
        std::max(1.0, std::max(std::abs(representative), std::abs(value)));

    if (std::abs(representative - value) > 1.0e-9 * scale) {

      return trimmed;
    }
  }

  if (haveRepresentative) {

    return formatTooltipScalar(representative);
  }

  return "0";
}

String componentTooltip(TopologicalPowerComp &component) {

  struct TooltipEntry {
    String name;
    String value;
  };

  std::vector<TooltipEntry> setpointEntries;

  std::vector<TooltipEntry> parameterEntries;

  std::vector<TooltipEntry> fallbackEntries;

  for (const auto &attribute : component.attributes()) {

    if (isRuntimeTooltipAttribute(attribute.first)) {

      continue;
    }

    TooltipEntry entry{
        attribute.first,
        compactTooltipValue(attribute.first, attribute.second->toString())};

    if (entry.value.empty()) {
      continue;
    }

    fallbackEntries.push_back(entry);

    if (isSetpointOrRatingAttribute(attribute.first)) {

      setpointEntries.push_back(entry);

      continue;
    }

    if (isLikelyConfiguredParameter(attribute.first)) {

      parameterEntries.push_back(std::move(entry));
    }
  }

  std::ostringstream tooltip;

  tooltip << component.name() << '\n' << component.type();

  if (!setpointEntries.empty()) {

    tooltip << "\nSetpoints / Ratings:";

    for (const auto &entry : setpointEntries) {

      tooltip << '\n' << entry.name << ": " << entry.value;
    }
  }

  if (!parameterEntries.empty()) {

    tooltip << "\nParameters:";

    for (const auto &entry : parameterEntries) {

      tooltip << '\n' << entry.name << ": " << entry.value;
    }
  }

  if (setpointEntries.empty() && parameterEntries.empty()) {

    tooltip << "\nParameters:";

    if (fallbackEntries.empty()) {

      tooltip << "\n(no configured parameters exposed as attributes)";

    } else {

      for (const auto &entry : fallbackEntries) {

        tooltip << '\n' << entry.name << ": " << entry.value;
      }
    }
  }

  return tooltip.str();
}

struct Dimensions {
  Real width = 0.40;
  Real height = 0.40;
};

Bool hasSymbolMapping(TopologicalPowerComp &component) {

  const String type = component.type();

  return contains(type, "GFM") || contains(type, "GFL") ||
         contains(type, "SynchronGenerator") || contains(type, "Transformer") ||
         contains(type, "PiLine") || contains(type, "Capacitor") ||
         contains(type, "Inductor") || contains(type, "Reactor") ||
         contains(type, "Resistor") || contains(type, "Switch");
}

Dimensions baseDimensions(TopologicalPowerComp &component) {

  const String type = component.type();

  if (contains(type, "SynchronGenerator"))
    return {0.40, 0.40};

  if (contains(type, "Transformer"))
    return {0.50, 0.36};

  if (contains(type, "PiLine"))
    return {0.58, 0.22};

  if (contains(type, "Resistor"))
    return {0.46, 0.22};

  if (contains(type, "Inductor") || contains(type, "Reactor"))
    return {0.46, 0.22};

  if (contains(type, "Capacitor"))
    return {0.38, 0.28};

  return {TopologyRenderDetail::RenderStyle::componentWidth,
          TopologyRenderDetail::RenderStyle::componentHeight};
}

#ifdef WITH_GRAPHVIZ

void configureGraph(Graph::Graph &graph,
                    const SystemTopologyRenderer::Options &options) {

  graph.set("bgcolor", "white");
  graph.set("outputorder", "edgesfirst");
  graph.set("splines", "ortho");
  graph.set("overlap", "false");
  graph.set("concentrate", "false");

  if (options.layout == SystemTopologyRenderer::Layout::LeftToRight) {

    graph.set("rankdir", "LR");

  } else {
    graph.set("rankdir", "TB");
  }

  Real nodeSeparation =
      std::max<Real>(TopologyRenderDetail::RenderStyle::nodeSeparation, 0.34);

  Real rankSeparation =
      std::max<Real>(TopologyRenderDetail::RenderStyle::rankSeparation, 0.34);

  if (options.layout == SystemTopologyRenderer::Layout::LeftToRight) {

    nodeSeparation = std::max<Real>(nodeSeparation, 0.48);

    rankSeparation = std::max<Real>(rankSeparation, 0.72);

  } else if (options.layout == SystemTopologyRenderer::Layout::TopToBottom) {

    nodeSeparation = std::max<Real>(nodeSeparation, 0.40);

    rankSeparation = std::max<Real>(rankSeparation, 0.44);

  } else {
    // Quadratic is only the Graphviz fallback for non-circular grids.
    graph.set("ratio", "compress");
    graph.set("pack", "true");
    graph.set("packmode", "graph");
  }

  graph.set("nodesep", std::to_string(nodeSeparation));

  graph.set("ranksep", std::to_string(rankSeparation));

  graph.set("esep", "+5");
  graph.set("sep", "+7");
  graph.set("margin", "0.06");
  graph.set("pad", "0.07");
}

void configureBus(Graph::Node &graphNode, const TopologicalNode::Ptr &node,
                  const SystemTopologyRenderer::Options &options) {

  const Bool horizontal =
      options.layout != SystemTopologyRenderer::Layout::LeftToRight;

  graphNode.set("shape", "box");
  graphNode.set("fixedsize", "true");

  graphNode.set("width", horizontal ? "0.52" : "0.075");

  graphNode.set("height", horizontal ? "0.075" : "0.52");

  graphNode.set("style", "filled");
  graphNode.set("color", "transparent");
  graphNode.set("fillcolor", "transparent");
  graphNode.set("penwidth", "0");
  graphNode.set("label", "");

  if (TopologyRenderDetail::RenderStyle::showNames) {
    graphNode.set("xlabel",
                  "<FONT FACE=\"" + String(kFont) + "\" POINT-SIZE=\"" +
                      formatReal(kBusLabelSize) + "\">" +
                      escapeLabel(node->name()) + "</FONT>",
                  true);
  }

  // Same convention as SystemTopology::topologyGraph(): bus hover identifies
  // the underlying DPsim node. Component hovers below expose all attributes.
  graphNode.set("tooltip", node->uid(), true);
}

void configureComponent(Graph::Node &graphNode,
                        TopologicalPowerComp::Ptr component,
                        const SystemTopologyRenderer::Options &options) {

  auto dimensions = baseDimensions(*component);

  const Bool verticalInitial =
      options.layout != SystemTopologyRenderer::Layout::LeftToRight;

  if (verticalInitial)
    std::swap(dimensions.width, dimensions.height);

  if (TopologyRenderDetail::RenderStyle::useSymbols &&
      hasSymbolMapping(*component)) {

    graphNode.set("shape", "box");
    graphNode.set("fixedsize", "true");

    graphNode.set("width", std::to_string(dimensions.width));

    graphNode.set("height", std::to_string(dimensions.height));

    graphNode.set("style", "filled");
    graphNode.set("color", "transparent");
    graphNode.set("fillcolor", "transparent");
    graphNode.set("penwidth", "0");
    graphNode.set("label", "");

  } else {
    // IMPORTANT:
    // Keep every generic/unmapped component as a plain Graphviz BOX.
    //
    // The topology renderer extracts and updates component geometry from the
    // polygon emitted by Graphviz. A rounded box may be emitted as an SVG
    // <path> instead of a <polygon>; such a component then disappears from the
    // geometry/orientation maps. This affected components such as Switch and
    // caused a later unordered_map::at() exception.
    graphNode.set("shape", "box");
    graphNode.set("fixedsize", "true");

    graphNode.set("width", std::to_string(dimensions.width));

    graphNode.set("height", std::to_string(dimensions.height));

    graphNode.set("style", "filled");
    graphNode.set("color", "black");
    graphNode.set("fillcolor", "white");
    graphNode.set("penwidth", "0.5");
    graphNode.set("fontcolor", "black");
    graphNode.set("fontname", kFont);
    graphNode.set("fontsize", "5.5");
    graphNode.set("margin", "0.03");

    graphNode.set("label", escapeLabel(component->type()));
  }

  if (TopologyRenderDetail::RenderStyle::showNames ||
      TopologyRenderDetail::RenderStyle::showTypes) {

    std::ostringstream label;

    label << "<FONT FACE=\"" << kFont << "\" POINT-SIZE=\""
          << kComponentLabelSize << "\">";

    if (TopologyRenderDetail::RenderStyle::showNames)
      label << escapeLabel(component->name());

    if (TopologyRenderDetail::RenderStyle::showNames &&
        TopologyRenderDetail::RenderStyle::showTypes)
      label << "<BR/>";

    if (TopologyRenderDetail::RenderStyle::showTypes) {
      label << "<FONT POINT-SIZE=\"" << kComponentTypeSize << "\">"
            << escapeLabel(component->type()) << "</FONT>";
    }

    label << "</FONT>";

    graphNode.set("xlabel", label.str(), true);
  }

  graphNode.set("tooltip", componentTooltip(*component), true);
}

std::unique_ptr<Graph::Graph>
buildGraph(const SystemTopology &system,
           const SystemTopologyRenderer::Options &options) {

  auto graph =
      std::make_unique<Graph::Graph>("dpsim_topology", Graph::Type::directed);

  configureGraph(*graph, options);

  std::unordered_map<const TopologicalNode *, Graph::Node *> graphNodes;

  for (const auto &node : system.mNodes) {
    if (!node || node->isGround())
      continue;

    auto *graphNode = graph->addNode("bus_" + node->uid());

    configureBus(*graphNode, node, options);

    graphNodes[node.get()] = graphNode;
  }

  UInt edgeCounter = 0;

  for (const auto &identifiedObject : system.mComponents) {
    if (!identifiedObject)
      continue;

    auto component =
        std::dynamic_pointer_cast<TopologicalPowerComp>(identifiedObject);

    if (!component)
      continue;

    const auto terminals = component->topologicalTerminals();

    // Optional collapsed line representation.
    if (TopologyRenderDetail::RenderStyle::collapseLines &&
        TopologyRenderDetail::isLineLike(*component) && terminals.size() == 2) {

      const auto node0 =
          terminals[0] ? terminals[0]->topologicalNodes() : nullptr;

      const auto node1 =
          terminals[1] ? terminals[1]->topologicalNodes() : nullptr;

      if (!node0 || !node1 || node0->isGround() || node1->isGround())
        continue;

      auto *edge = graph->addEdge("line_" + std::to_string(edgeCounter++),
                                  graphNodes.at(node0.get()),
                                  graphNodes.at(node1.get()));

      edge->set("arrowhead", "none");
      edge->set("color", "black");
      edge->set("penwidth", "0.25");
      edge->set("tooltip", componentTooltip(*component), true);

      if (TopologyRenderDetail::RenderStyle::showNames) {
        edge->set("xlabel",
                  "<FONT FACE=\"" + String(kFont) + "\" POINT-SIZE=\"" +
                      formatReal(kLineLabelSize) + "\">" +
                      escapeLabel(component->name()) + "</FONT>",
                  true);
      }

      continue;
    }

    auto *componentNode = graph->addNode("component_" + component->uid());

    configureComponent(*componentNode, component, options);

    if (terminals.size() == 1) {
      if (!terminals[0])
        continue;

      const auto node = terminals[0]->topologicalNodes();

      if (!node || node->isGround())
        continue;

      Graph::Edge *edge = nullptr;

      // Keep sources/generators on the "outside" of Graphviz's initial rank.
      if (contains(component->type(), "SynchronGenerator")) {

        edge = graph->addEdge("connection_" + std::to_string(edgeCounter++),
                              componentNode, graphNodes.at(node.get()));

      } else {
        edge = graph->addEdge("connection_" + std::to_string(edgeCounter++),
                              graphNodes.at(node.get()), componentNode);
      }

      edge->set("arrowhead", "none");
      edge->set("color", "black");
      edge->set("penwidth", "0.25");
      continue;
    }

    if (terminals.size() == 2) {
      const auto node0 =
          terminals[0] ? terminals[0]->topologicalNodes() : nullptr;

      const auto node1 =
          terminals[1] ? terminals[1]->topologicalNodes() : nullptr;

      const Bool ground0 = !node0 || node0->isGround();

      const Bool ground1 = !node1 || node1->isGround();

      if (!ground0 && !ground1) {
        auto *edge0 =
            graph->addEdge("connection_" + std::to_string(edgeCounter++),
                           graphNodes.at(node0.get()), componentNode);

        edge0->set("arrowhead", "none");
        edge0->set("color", "black");
        edge0->set("penwidth", "0.25");

        auto *edge1 =
            graph->addEdge("connection_" + std::to_string(edgeCounter++),
                           componentNode, graphNodes.at(node1.get()));

        edge1->set("arrowhead", "none");
        edge1->set("color", "black");
        edge1->set("penwidth", "0.25");
        continue;
      }

      const auto networkNode = ground0 ? node1 : node0;

      if (!networkNode || networkNode->isGround())
        continue;

      auto *edge =
          graph->addEdge("shunt_" + std::to_string(edgeCounter++),
                         graphNodes.at(networkNode.get()), componentNode);

      edge->set("arrowhead", "none");
      edge->set("color", "black");
      edge->set("penwidth", "0.25");

      // Shunts never stretch the electrical backbone.
      edge->set("constraint", "false");
      continue;
    }

    // Generic multi-terminal fallback.
    for (const auto &terminal : terminals) {
      if (!terminal)
        continue;

      const auto node = terminal->topologicalNodes();

      if (!node || node->isGround())
        continue;

      auto *edge = graph->addEdge("connection_" + std::to_string(edgeCounter++),
                                  graphNodes.at(node.get()), componentNode);

      edge->set("arrowhead", "none");
      edge->set("color", "black");
      edge->set("penwidth", "0.25");
    }
  }

  return graph;
}

#endif // WITH_GRAPHVIZ

} // namespace

SystemTopologyRenderer::SystemTopologyRenderer(
    const SystemTopology &system, std::filesystem::path symbolDirectory)
    : SystemTopologyRenderer(system, std::move(symbolDirectory), Options{}) {}

SystemTopologyRenderer::SystemTopologyRenderer(
    const SystemTopology &system, std::filesystem::path symbolDirectory,
    Options options)
    : mSystem(system), mSymbolDirectory(std::move(symbolDirectory)),
      mOptions(options) {

  if (TopologyRenderDetail::RenderStyle::useSymbols &&
      !std::filesystem::exists(mSymbolDirectory)) {

    throw std::runtime_error(
        "Topology symbol directory does not exist: " +
        std::filesystem::absolute(mSymbolDirectory).string());
  }
}

void SystemTopologyRenderer::renderSvg(
    const std::filesystem::path &filename) const {

#ifndef WITH_GRAPHVIZ

  throw std::runtime_error("SystemTopologyRenderer requires WITH_GRAPHVIZ");

#else

  if (!filename.parent_path().empty()) {
    std::filesystem::create_directories(filename.parent_path());
  }

  std::ofstream stream(filename);

  if (!stream.is_open()) {
    throw std::runtime_error("Cannot open topology SVG output file: " +
                             filename.string());
  }

  stream << renderSvg();

#endif
}

String SystemTopologyRenderer::renderSvg() const {

#ifndef WITH_GRAPHVIZ

  throw std::runtime_error("SystemTopologyRenderer requires WITH_GRAPHVIZ");

#else

  if (TopologyRenderDetail::RenderStyle::useSymbols &&
      !std::filesystem::exists(mSymbolDirectory)) {

    throw std::runtime_error(
        "Topology symbol directory does not exist: " +
        std::filesystem::absolute(mSymbolDirectory).string());
  }

  auto graph = buildGraph(mSystem, mOptions);

  std::stringstream graphStream;

  const auto graphvizStart = std::chrono::steady_clock::now();

  graph->render(graphStream, "dot", "svg");

  const auto graphvizEnd = std::chrono::steady_clock::now();

  String svg = graphStream.str();

  // ---------------------------------------------------------------------------
  // Layout algorithm
  //
  // 1. Detect simple closed electrical ring.
  //    -> arrange the ring first;
  //    -> place every external branch outside the ring.
  //
  // 2. If no ring:
  //    LeftToRight / TopToBottom
  //    -> find the longest bus-to-bus connection;
  //    -> use it as the publication backbone;
  //    -> everything connected to a middle backbone bus becomes a leaf.
  //
  // 3. Quadratic without a simple ring:
  //    -> retain Graphviz's compact fallback.
  // ---------------------------------------------------------------------------

  const auto layoutStart = std::chrono::steady_clock::now();

  const auto layoutResult = TopologyRenderDetail::applyTopologyLayout(
      svg, mSystem, mOptions.layout,
      TopologyRenderDetail::RenderStyle::collapseLines);

  const auto layoutEnd = std::chrono::steady_clock::now();

  String result = TopologyRenderDetail::postProcessSvg(
      std::move(svg), mSystem, mSymbolDirectory, mOptions, layoutResult);

  const auto graphvizMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(graphvizEnd -
                                                            graphvizStart)
          .count();

  const auto layoutMilliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(layoutEnd -
                                                            layoutStart)
          .count();

  SPDLOG_INFO(
      "Topology renderer stages: Graphviz={} ms, topology-layout={} ms.",
      graphvizMilliseconds, layoutMilliseconds);

  return result;

#endif
}

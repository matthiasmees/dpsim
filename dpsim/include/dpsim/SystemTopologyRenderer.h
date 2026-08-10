// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/SystemTopology.h>

#include <filesystem>

namespace CPS {

class SystemTopologyRenderer {
public:
  enum class Layout { LeftToRight, TopToBottom, Quadratic };

  // Public configuration is intentionally restricted to the layout choice.
  // All visual styling, spacing, collision handling and routing parameters are
  // renderer-internal defaults.
  struct Options {
    Layout layout = Layout::LeftToRight;
  };

  SystemTopologyRenderer(const SystemTopology &system,
                         std::filesystem::path symbolDirectory);

  SystemTopologyRenderer(const SystemTopology &system,
                         std::filesystem::path symbolDirectory,
                         Options options);

  void renderSvg(const std::filesystem::path &filename) const;

  String renderSvg() const;

private:
  const SystemTopology &mSystem;
  std::filesystem::path mSymbolDirectory;
  Options mOptions;
};

} // namespace CPS

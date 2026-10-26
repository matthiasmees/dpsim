/* A data logger for real-time simulation data logging.
 *
 * Author: Niklas Eiling <niklas.eiling@eonerc.rwth-aachen.de>
 * SPDX-FileCopyrightText: 2024 Niklas Eiling <niklas.eiling@eonerc.rwth-aachen.de>
 * SPDX-License-Identifier: Apache-2.0
 */

#include <algorithm>
#include <iomanip>
#include <memory>

#include <dpsim/RealTimeDataLogger.h>

#include <dpsim-models/Attribute.h>
#include <dpsim-models/Logger.h>

using namespace DPsim;

RealTimeDataLogger::RealTimeDataLogger(std::filesystem::path &filename,
                                       size_t rowNumber)
    : DataLoggerInterface(), mFilename(filename), mRowNumber(rowNumber),
      mCurrentRow(0), mCurrentAttribute(0), mRowsWritten(0), mAttributeData() {}

RealTimeDataLogger::RealTimeDataLogger(std::filesystem::path &filename,
                                       Real finalTime, Real timeStep)
    : DataLoggerInterface(), mFilename(filename),
      // Include the initial t=0 row and reserve one more row for an endpoint
      // that can be reached exactly by the simulation's floating-point loop.
      mRowNumber(static_cast<size_t>(finalTime / timeStep + 0.5) + 1),
      mCurrentRow(0), mCurrentAttribute(0), mRowsWritten(0), mAttributeData() {}

void RealTimeDataLogger::start() {
  double mb_size =
      static_cast<double>(mRowNumber) * (mAttributes.size() + 1) * sizeof(Real);
  auto log = CPS::Logger::get("RealTimeDataLogger", CPS::Logger::Level::off,
                              CPS::Logger::Level::info);
  SPDLOG_LOGGER_INFO(
      log,
      "Preallocating memory for real-time data logger: {} rows for {} "
      "attributes ({} MB)",
      mRowNumber, mAttributes.size(), mb_size / (1024 * 1024));
  // We are doing real time so preallocate everything
  mAttributeData.resize(mRowNumber);
  for (auto &it : mAttributeData) {
    // We have to add one to the size because we also log the time
    it.resize(mAttributes.size() + 1);
  }
}

void RealTimeDataLogger::stop() {
  auto log = CPS::Logger::get("RealTimeDataLogger", CPS::Logger::Level::off,
                              CPS::Logger::Level::info);
  SPDLOG_LOGGER_INFO(
      log, "Stopping real-time data logger. Writing memory to file {}",
      mFilename.string());

  const auto parent = mFilename.parent_path();
  if (!parent.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(parent, ec);
    if (ec) {
      throw std::runtime_error("Cannot create log directory '" +
                               parent.string() + "': " + ec.message());
    }
  }

  auto mLogFile =
      std::ofstream(mFilename, std::ios_base::out | std::ios_base::trunc);
  if (!mLogFile.is_open()) {
    throw std::runtime_error("Cannot open log file " + mFilename.string());
  }

  mLogFile << std::right << std::setw(14) << "time";
  for (auto it : mAttributes)
    mLogFile << ", " << std::right << std::setw(13) << it.first;
  mLogFile << '\n';

  for (size_t rowIdx = 0; rowIdx < mRowsWritten; ++rowIdx) {
    const auto &row = mAttributeData[rowIdx];
    mLogFile << std::scientific << std::right << std::setw(14) << row[0];
    for (size_t i = 1; i < row.size(); ++i)
      mLogFile << ", " << std::fixed << std::setprecision(6) << std::right
               << std::setw(13) << row[i];
    mLogFile << '\n';
  }
  mLogFile.close();
  SPDLOG_LOGGER_INFO(log, "Finished writing real-time data log to file {}",
                     mFilename.string());
}

void RealTimeDataLogger::log(Real time, Int timeStepCount) {
  // Simulation::start logs the t=0 state with count 0; the first solved step
  // also has count 0, but occurs at t=dt. Keep both rows by offsetting all
  // positive-time samples by one.
  if (timeStepCount < 0) {
    throw std::runtime_error(
        "RealTimeDataLogger: timeStepCount out of bounds. Please verify the "
        "logger was initialized correctly.");
  }
  const size_t row = static_cast<size_t>(timeStepCount) + (time > 0 ? 1 : 0);
  mCurrentRow = row;
  if (row >= mRowNumber) {
    throw std::runtime_error(
        "RealTimeDataLogger: timeStepCount out of bounds. Please verify the "
        "logger was initialized correctly.");
  }
  if (mAttributeData.size() != mRowNumber ||
      mAttributeData[row].size() != mAttributes.size() + 1) {
    throw std::runtime_error(
        "RealTimeDataLogger: Attribute data size mismatch");
  }
  mAttributeData[row][0] = time;
  mRowsWritten = std::max(mRowsWritten, row + 1);
  mCurrentAttribute = 1;

  for (auto &it : mAttributes) {
    auto base = it.second.getPtr(); // std::shared_ptr<CPS::AttributeBase>

    if (it.second->getType() == typeid(Real)) {
      auto attr = std::dynamic_pointer_cast<CPS::Attribute<Real>>(base);
      if (!attr) { /* skip */
        continue;
      }
      mAttributeData[row][mCurrentAttribute++] = attr->get();
    } else if (it.second->getType() == typeid(Int)) {
      auto attr = std::dynamic_pointer_cast<CPS::Attribute<Int>>(base);
      if (!attr) { /* skip */
        continue;
      }
      mAttributeData[row][mCurrentAttribute++] = attr->get();
    }
  }
}

void RealTimeDataLogger::Step::execute(Real time, Int timeStepCount) {
  mLogger.log(time, timeStepCount);
}

CPS::Task::Ptr RealTimeDataLogger::getTask() {
  return std::make_shared<RealTimeDataLogger::Step>(*this);
}

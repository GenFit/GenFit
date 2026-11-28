/* Copyright 2015, Ludwig-Maximilians-Universität München
   Authors: Tobias Schlüter

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

/** @brief Defines for I/O streams used for error and debug printing.
 */

#include <mutex>
#include <ostream>

namespace genfit {

/** @brief Output stream whose insertions are serialized between fit threads.
 *
 * Every operator<< holds a mutex for the duration of the insertion, so several
 * threads may write through the same stream without interleaving or racing on the
 * stream state. Direct std::ostream access is not provided on purpose; swap the
 * destination with rdbuf() instead.
 */
class SyncStream {
 public:
  explicit SyncStream(std::streambuf* buf) : stream_(buf) {}

  template <typename T>
  SyncStream& operator<<(const T& value) {
    std::lock_guard<std::mutex> lock(mutex_);
    stream_ << value;
    return *this;
  }

  //! Manipulators, e.g. std::endl (a function template, not deducible above).
  SyncStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
    std::lock_guard<std::mutex> lock(mutex_);
    manip(stream_);
    return *this;
  }
  SyncStream& operator<<(std::ios& (*manip)(std::ios&)) {
    std::lock_guard<std::mutex> lock(mutex_);
    manip(stream_);
    return *this;
  }
  SyncStream& operator<<(std::ios_base& (*manip)(std::ios_base&)) {
    std::lock_guard<std::mutex> lock(mutex_);
    manip(stream_);
    return *this;
  }

  //! Same interface as std::ostream::rdbuf(): get or swap the destination buffer.
  std::streambuf* rdbuf() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stream_.rdbuf();
  }
  std::streambuf* rdbuf(std::streambuf* buf) {
    std::lock_guard<std::mutex> lock(mutex_);
    return stream_.rdbuf(buf);
  }

 private:
  SyncStream(const SyncStream&) = delete;
  SyncStream& operator=(const SyncStream&) = delete;

  mutable std::mutex mutex_;
  std::ostream stream_;
};

/** Default stream for debug output.  Defaults to std::cout.
   Override destination with debugOut.rdbuf(newStream.rdbuf()).  */
extern SyncStream debugOut;
/** Default stream for error output.  Defaults to std::cerr.
    Override destination with errorOut.rdbuf(newStream.rdbuf()).  */
extern SyncStream errorOut;
/** Default stream for output of Print calls.  Defaults to std::cout.
   Override destination with printOut.rdbuf(newStream.rdbuf()).  */
extern SyncStream printOut;

}

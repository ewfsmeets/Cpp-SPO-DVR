// gfortran sequential unformatted records (length markers, subrecords).
#include "io/fortran_unformatted.hpp"

#include <algorithm>
#include <cstdlib>

namespace spodvr {

void FortranRecordWriter::write_bytes(std::span<const std::byte> data) {
  std::size_t offset = 0;
  bool first = true;
  do {
    const auto length = static_cast<std::int64_t>(std::min<std::size_t>(data.size() - offset, max_subrecord_));
    const bool more = offset + length < data.size();
    const auto marker = static_cast<std::int32_t>(length);
    file_.write<std::int32_t>(more ? -marker : marker);  // head: negative if continued later
    file_.write_array(data.subspan(offset, length));
    file_.write<std::int32_t>(first ? marker : -marker);  // tail: negative if a continuation
    offset = offset + length;
    first = false;
  } while (offset < data.size());
}

void FortranRecordParser::read_bytes(std::span<std::byte> out) {
  if (position_ + out.size() > bytes_.size())
    throw CheckpointError(description_ + ": record too short (" + std::to_string(bytes_.size()) + " bytes)");
  std::copy_n(bytes_.begin() + static_cast<std::ptrdiff_t>(position_), out.size(), out.begin());
  position_ = position_ + out.size();
}

void FortranRecordParser::expect_end() const {
  if (position_ != bytes_.size())
    throw CheckpointError(description_ + ": record has " + std::to_string(bytes_.size()) + " bytes, expected " +
                          std::to_string(position_));
}

FortranRecordParser FortranRecordReader::read(const std::string& what) {
  const std::string description = what + " in '" + file_.path().string() + "'";
  std::vector<std::byte> bytes;
  bool more = true;
  while (more) {
    const auto head = file_.read<std::int32_t>();
    const std::size_t length = std::abs(static_cast<std::int64_t>(head));
    const std::size_t offset = bytes.size();
    bytes.resize(offset + length);
    file_.read_array(std::span(bytes).subspan(offset, length));
    const auto tail = file_.read<std::int32_t>();
    if (std::abs(static_cast<std::int64_t>(tail)) != static_cast<std::int64_t>(length))
      throw CheckpointError(description + ": inconsistent record markers (not a gfortran unformatted file?)");
    more = head < 0;
  }
  return FortranRecordParser(std::move(bytes), description);
}

}  // namespace spodvr

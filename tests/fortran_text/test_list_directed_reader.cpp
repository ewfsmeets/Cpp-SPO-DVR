// Tests: list-directed reading of text files, the subset the CRP data-set files use
// (src/fortran_text/list_directed_reader.*).
#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "fortran_text/list_directed_reader.hpp"
#include "temporary_directory.hpp"

using fortran_text::ListDirectedReader;
using fortran_text::parse_fortran_real;
using fortran_text::ReadError;
using test_support::TemporaryDirectory;

namespace {

ListDirectedReader reader_for(const TemporaryDirectory& dir, std::string_view text) {
  return ListDirectedReader(dir.write("input.dat", text));
}

}  // namespace

TEST(FortranReal, AcceptsTheFortranForms) {
  EXPECT_EQ(parse_fortran_real("1"), 1.0);
  EXPECT_EQ(parse_fortran_real("1."), 1.0);
  EXPECT_EQ(parse_fortran_real(".5"), 0.5);
  EXPECT_EQ(parse_fortran_real("+2.5"), 2.5);
  EXPECT_EQ(parse_fortran_real("-0.25"), -0.25);
  EXPECT_EQ(parse_fortran_real("1.5e3"), 1500.0);
  EXPECT_EQ(parse_fortran_real("1.5E+3"), 1500.0);
  EXPECT_EQ(parse_fortran_real("-1.0d-3"), -1.0e-3);
  EXPECT_EQ(parse_fortran_real("2D2"), 200.0);
  EXPECT_EQ(parse_fortran_real("0.529177249"), 0.529177249);
}

TEST(FortranReal, IsCorrectlyRounded) {
  // More digits than a double holds: the nearest double (as gfortran's strtod).
  EXPECT_EQ(parse_fortran_real("3.1415926535897932384626433832795028841971"), 3.141592653589793);
  EXPECT_EQ(parse_fortran_real("0.1"), 0.1);
  EXPECT_EQ(parse_fortran_real("1.0000000000000001110223024625156540423631668090820312500000000001"),
            std::nextafter(1.0, 2.0));
}

TEST(FortranReal, RejectsOtherText) {
  for (const char* text : {"", "+", ".", "e5", "1e", "1.0+3", "1q0", "abc", "1.0x", "inf", "nan", "1,0", "--1"})
    EXPECT_FALSE(parse_fortran_real(text).has_value()) << text;
}

TEST(ListDirectedReader, OneStatementPerRecord) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "C6v  extra words\n1 3\n4.8\n");
  EXPECT_EQ(in.read(1).text(0), "C6v");  // the rest of the record is skipped
  const auto units = in.read(2);
  EXPECT_EQ(units.integer(0), 1);
  EXPECT_EQ(units.integer(1), 3);
  EXPECT_EQ(in.read(1).real(0), 4.8);
  EXPECT_FALSE(in.try_read(1).has_value());
  EXPECT_THROW(in.read(1), ReadError);
}

TEST(ListDirectedReader, StatementContinuesOnTheNextRecords) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "1.0\n2.0, 3.0\n\n4.0 5.0 6.0\n7.0\n");
  const auto s = in.read(4);
  ASSERT_EQ(s.size(), 4);
  EXPECT_EQ(s.real(3), 4.0);
  EXPECT_EQ(in.read(1).real(0), 7.0);  // 5.0 and 6.0 were skipped with the rest of the record
}

TEST(ListDirectedReader, SeparatorsAndQuotes) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "\t1,2 ,  3\r\n'a b.dat' \"it''s\" 'x''y'\n");
  const auto numbers = in.read(3);
  EXPECT_EQ(numbers.integer(2), 3);
  const auto names = in.read(3);
  EXPECT_EQ(names.text(0), "a b.dat");
  EXPECT_EQ(names.text(1), "it''s");  // doubled quotes only count inside the same kind of quotes
  EXPECT_EQ(names.text(2), "x'y");
}

TEST(ListDirectedReader, LogicalValues) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "18 T\n18 .FALSE.\n1 .true.\n1 f\n1 yes\n");
  EXPECT_TRUE(in.read(2).logical(1));
  EXPECT_FALSE(in.read(2).logical(1));
  EXPECT_TRUE(in.read(2).logical(1));
  EXPECT_FALSE(in.read(2).logical(1));
  EXPECT_THROW((void)in.read(2).logical(1), ReadError);
}

TEST(ListDirectedReader, EmptyStatementSkipsARecord) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "1\nskipped line\n2.5 3.5\n");
  EXPECT_EQ(in.read(1).integer(0), 1);
  EXPECT_EQ(in.read(0).size(), 0);  // e.g. NumLayers - 1 = 0 interlayer spacings
  EXPECT_EQ(in.read(2).real(1), 3.5);
}

TEST(ListDirectedReader, OptionalStatementAtTheEnd) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "1\n0.5\n");
  EXPECT_EQ(in.read(1).integer(0), 1);
  EXPECT_FALSE(in.try_read(2).has_value());  // one value only, then the end of the file (END=)
}

TEST(ListDirectedReader, LinesAsTheyAre) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "  short info, with / and *\r\nlong\n");
  EXPECT_EQ(in.read_line(), "  short info, with / and *");
  EXPECT_EQ(in.read_line(), "long");
  EXPECT_FALSE(in.read_line().has_value());
}

TEST(ListDirectedReader, ErrorsNameFileAndLine) {
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "1 2\n3*1.0\n");
  (void)in.read(2);
  try {
    (void)in.read(3);
    FAIL() << "repeat count accepted";
  } catch (const ReadError& e) {
    const std::string message = e.what();
    EXPECT_NE(message.find("input.dat:2:"), std::string::npos) << message;
    EXPECT_NE(message.find("repeat counts"), std::string::npos) << message;
  }
}

TEST(ListDirectedReader, UnsupportedInputIsRejected) {
  TemporaryDirectory dir("reader");
  EXPECT_THROW((void)reader_for(dir, "1,,2\n").read(3), ReadError);      // null value
  EXPECT_THROW((void)reader_for(dir, ",1 2\n").read(2), ReadError);      // leading null value
  EXPECT_THROW((void)reader_for(dir, "1 / 2\n").read(2), ReadError);     // slash
  EXPECT_THROW((void)reader_for(dir, "'open\n").read(1), ReadError);     // unterminated string
  EXPECT_THROW((void)reader_for(dir, "1.5\n").read(1).integer(0), ReadError);  // not an integer
  EXPECT_THROW((void)reader_for(dir, "x\n").read(1).real(0), ReadError);
  EXPECT_THROW(ListDirectedReader(dir.path() / "missing.dat"), ReadError);
}

TEST(ListDirectedReader, TextAfterTheLastValueIsIgnored) {
  // Fortran stops reading a record once it has all values: the rest may contain anything.
  TemporaryDirectory dir("reader");
  auto in = reader_for(dir, "1.0 2.0 / comment, 3*4 ,,\n");
  EXPECT_EQ(in.read(2).real(1), 2.0);
}

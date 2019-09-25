// This file is part of posrecon.
// 
// posrecon is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// posrecon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with posrecon.  If not, see <http://www.gnu.org/licenses/>.
//
// Written by SMK <chairmanK@gmail.com>. Last revised 2009-02-01.
//
// VideoReader.h defines a VideoReader class for posrecon. This class handles 
// the following NSpike data files:
//    * MPEG-2 file that consists entirely of I-frames. Non-compliant MPEG-2 
//      files will not be read correctly! 
//    * mpeg index file, which contains filepos_t byte positions into the mpeg 
//      file, allowing for fast seek-to-frame functionality.
//    * timestamp file, which contains uint32_t timestamps that match the
//      frames of the mpeg video.

#include <stdint.h>
extern "C" {
  #include "mpeg2dec/mpeg2.h"
  #include "mpeg2dec/mpeg2convert.h"
}
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#ifndef VIDEOREADER_H_
#define VIDEOREADER_H_

typedef uint32_t framecount_t;
typedef int64_t frameoff_t;
typedef uint64_t filepos_t;

// Opens mpeg file, mpeg index file, timestamps file. Reads image 
// frames and matching timestamps. 
class VideoReader {
 public:
  VideoReader();
  ~VideoReader();
  void open(const std::string& mpeg_stream_filename, 
      const std::string& mpeg_index_filename, 
      const std::string& timestamps_filename); 
  void close();
  const bool& ready() const;
  const unsigned int& width() const;
  const unsigned int& height() const;
  const unsigned int& aspect_ratio_numerator() const;
  const unsigned int& aspect_ratio_denominator() const;
  const mpeg2convert_rgb_order_t& rgb_order() const;
  const unsigned int& bits_per_pixel() const;
  const framecount_t& number_of_frames() const;
  const framecount_t& current_index() const;
  const uint32_t& current_timestamp() const;
  /* returns all timestamps in sequential order as a C++ STL vector */
  const std::vector<uint32_t>& all_timestamps() const;
  /* returns a pointer to an array of uint8_t RGB values for the pixels in the
  last read video frame */
  const uint8_t* current_picture() const;
  /* seek_frame method is polymorphic: if given a stream offset and origin, then
  it goes to that frame (exactly) or returns false if no such frame exists. if
  given a timestamp argument, then it goes to the nearest frame or returns
  false if the specified timestamp is outside of the range of timestamps */
  bool seek_frame(frameoff_t frame_offset, std::ios::seekdir dir);
  bool seek_frame(uint32_t timestamp);

 private:
  VideoReader(const VideoReader&); // disable copy
  void operator=(const VideoReader&); // disable assignment
  bool open_mpeg_stream(const std::string& mpeg_stream_filename);
  bool open_mpeg_index(const std::string& mpeg_index_filename);
  bool open_timestamps(const std::string& timestamps_filename);
  mpeg2dec_t* mpeg_decoder_;
  const mpeg2_info_t* mpeg_info_; 
  bool ready_;
  framecount_t number_of_frames_;
  unsigned int width_;
  unsigned int height_;
  unsigned int aspect_ratio_numerator_;
  unsigned int aspect_ratio_denominator_;
  std::ifstream mpeg_stream_;
  std::vector<filepos_t> mpeg_index_;
  std::vector<uint32_t> timestamps_;
  uint32_t current_timestamp_;
  framecount_t current_index_;
  /* Before the mpeg decoder is ready, the framebuffer array does not exist.
  Hence we can't declare current_picture_ as a const reference. Instead, it is
  initialized as null pointer, which gets reseated once the decoded picture
  becomes available */
  const uint8_t* current_picture_;
};

#endif // VIDEOREADER_H_

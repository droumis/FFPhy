
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

#include "VideoReader.h"


/* set this constant large enough so that an entire frame fits into a single
 read buffer */
const size_t MPEG_BUFFER_SIZE = 230400;
/* PIXEL_CONVERSION_SYMBOL is a pointer to the desired mpeg conversion function;
 we want packed 8-bit RGB */
const mpeg2_convert_t* PIXEL_CONVERSION_SYMBOL = mpeg2convert_rgb24;
const mpeg2convert_rgb_order_t RGB_ORDER = MPEG2CONVERT_RGB;
const unsigned int BITS_PER_PIXEL = 8;


VideoReader::VideoReader() 
    : mpeg_decoder_(mpeg2_init()),
      mpeg_info_(0),
      mpeg_index_(),
      timestamps_(),
      ready_(false),
      number_of_frames_(0),
      width_(0),
      height_(0),
      current_index_(0),
      current_timestamp_(0),
      current_picture_(0),
      aspect_ratio_numerator_(0),
      aspect_ratio_denominator_(0),
      mpeg_stream_() {
}


void VideoReader::open(const std::string& mpeg_stream_filename, 
    const std::string& mpeg_index_filename, 
    const std::string& timestamps_filename) {
  if (ready_) {
    std::cerr << "Video reader object is already initialized\n";
    return;
  }
  /* Populate timestamps_ */
  if (!open_timestamps(timestamps_filename)) {
    std::cerr << "Failed to open timestamps\n";
    return;
  }
  /* Populate mpeg_index_ (vector of file stream offsets into the mpeg file) */
  if (!open_mpeg_index(mpeg_index_filename)) {
    std::cerr << "Failed to open mpeg index\n";
    return;
  }
  /* Check that mpeg_stream_ is open for reading and populate mpeg decoder
  object with initial data */
  if (!open_mpeg_stream(mpeg_stream_filename)) {
    std::cerr << "Failed to open mpeg stream\n";
    return;
  }
  /* Check that timestamps_ and mpeg_index_ match in size */
  if (mpeg_index_.size() == timestamps_.size()) {
    number_of_frames_ = static_cast<framecount_t>(mpeg_index_.size());
  } else {
    std::cerr << "Number of timestamps (" << timestamps_.size() << 
        ") does not match number of mpeg file positions(" << 
        mpeg_index_.size() << '\n';
  } 
  /* Seek to the first frame. As a side effect, this populates
  mpeg_info_->sequence) */
  if (!seek_frame(0, std::ifstream::beg)) {
    std::cerr << "Could not decode the first frame of the mpeg\n";
    return;
  }
  width_ = mpeg_info_->sequence->width;
  height_ = mpeg_info_->sequence->height;
  aspect_ratio_numerator_ = mpeg_info_->sequence->pixel_width;
  aspect_ratio_denominator_ = mpeg_info_->sequence->pixel_height;
  /* Seek to the last frame of the mpeg (to verify that the last element of
  mpeg_index_ is valid) */
  if (!seek_frame(-1, std::ifstream::end)) {
    std::cerr << "Could not decode the last frame of the mpeg\n";
    return;
  }
  /* Report information */
  std::cout << "successfully initialized.\npicture size is " << width_ << 
      "x" << height_ << " pixels\n" << "pixel aspect ratio is " << 
      aspect_ratio_numerator_ << ":" << aspect_ratio_denominator_ << 
      "\n" << "total number of frames is " << mpeg_index_.size() << "\n";
  ready_ = true;
}


void VideoReader::close() {
  mpeg_info_ = 0;
  number_of_frames_ = 0;
  mpeg_index_.resize(number_of_frames_);
  timestamps_.resize(number_of_frames_);
  width_ = 0;
  height_ = 0;
  current_index_ = 0;
  current_timestamp_ = 0;
  current_picture_ = 0;
  aspect_ratio_numerator_ = 0;
  aspect_ratio_denominator_ = 0;
  mpeg_stream_.close();
  ready_ = false;
}


VideoReader::~VideoReader() {
  close();
  mpeg2_close(mpeg_decoder_);
}


bool VideoReader::open_mpeg_stream(const std::string& mpeg_stream_filename) {
  // Guard against calling this function more than once  
  if (ready_) {
    std::cerr << "Object is already initialized\n";
    return false;
  }
  // Open mpeg stream
  mpeg_stream_.open(mpeg_stream_filename.c_str(), std::ifstream::binary);
  if (!mpeg_stream_.is_open()) {
    std::cerr << "mpeg file " << mpeg_stream_filename << "is not readable\n";
    return false;
  }
  // Check that the mpeg decoder was properly initialized
  if (mpeg_decoder_ == 0) {
    std::cerr << "Could not create an mpeg decoder instance\n";
    return false;
  }
  // Set mpeg2_info pointer
  mpeg_info_ = mpeg2_info(mpeg_decoder_);
  if (mpeg_info_ == 0) {
    std::cerr << "Could not read info from the mpeg decoder object\n";
    return false;
  }
  // Reset mpeg_decoder_ (we set the full_reset argument to 0 because we only
  // seek within the same mpeg sequence)
  mpeg2_reset(mpeg_decoder_, 0);
  /* std::cout << "reset mpeg decoder\n"; */
  
  /* std::cout << "initialized mpeg stream\n"; */
  return true;
}


bool VideoReader::open_mpeg_index(const std::string& mpeg_index_filename) {
  // Guard against calling this function more than once  
  if (ready_) {
    std::cerr << "Object is already initialized\n";
    return false;
  }
  // Open mpeg index file
  std::ifstream mpeg_index_ifstream(mpeg_index_filename.c_str());
  if (!mpeg_index_ifstream.is_open()) {
    std::cerr << "Could not open mpeg index file " << mpeg_index_filename << 
        '\n'; 
    return false;
  }
  // Scan lines of the mpeg index file for the EOH sentinel
  mpeg_index_ifstream.seekg(0);
  std::string line;
  std::string eoh_sentinel("%%ENDHEADER");
  std::streampos eoh_pos = -1;
  while (!mpeg_index_ifstream.eof()) {
    std::getline(mpeg_index_ifstream, line, '\n');
    // '%%ENDHEADER' must appear as its own line at EOH
    if ((line.find(eoh_sentinel) != std::string::npos) &&
        (eoh_sentinel.length() == line.length())) {
      eoh_pos = mpeg_index_ifstream.tellg();
      break;
    }
  }
  if (eoh_pos < 0) {
    std::cerr << "File " << mpeg_index_filename << 
        " is missing end-of-header line " << eoh_sentinel << '\n';
    return false;
  }
  // Re-open the mpeg index file as binary to read the data
  mpeg_index_ifstream.close();
  mpeg_index_ifstream.open(mpeg_index_filename.c_str(), std::ifstream::binary);
  // Estimate the size of the data block
  mpeg_index_ifstream.seekg(0, std::ifstream::end);
  std::streampos data_size = mpeg_index_ifstream.tellg() - eoh_pos;
  // Seek to the putative start of the binary data block after header
  mpeg_index_ifstream.seekg(eoh_pos);
  if (mpeg_index_ifstream.rdstate()) {
    std::cerr << "Can not seek to position " << eoh_pos << " in " << 
        mpeg_index_filename << '\n';
    return false;
  }
  // Resize vector to match the number of (complete) mpeg file positions
  mpeg_index_.resize(data_size/sizeof(filepos_t));
  // Read bytes and put them at &mpeg_index_ (is this dangerous?)
  if ((sizeof(filepos_t)*mpeg_index_.size() == data_size) & (data_size > 0)) {
    // (is there a more elegant way to do this with iterator?)
    mpeg_index_ifstream.read((char*) &mpeg_index_[0], data_size); 
  }
  mpeg_index_ifstream.close();
  // Check that mpeg_index_ are monotonic (remember that the NSpike mpeg file
  // ought to consist entirely of I-frames)
  std::vector<filepos_t>::iterator mpeg_index_iter = mpeg_index_.begin();
  while (mpeg_index_iter != mpeg_index_.end()) {
    if (*(mpeg_index_iter+1) - *(mpeg_index_iter) <= 0) {
      std::cerr << "File positions in mpeg index are out of order\n";
      return false;
    }
    mpeg_index_iter++;
  }
  /* std::cout << "initialized mpeg index\n"; */
  return true;
}


bool VideoReader::open_timestamps(const std::string& timestamps_filename) {
  // Guard against calling this function more than once  
  if (ready_) {
    std::cerr << "Object is already initialized\n";
    return false;
  }
  // Open timestamps file
  std::ifstream timestamps_ifstream(timestamps_filename.c_str());
  if (!timestamps_ifstream.is_open()) {
    std::cerr << "Could not open timestamps file " << timestamps_filename << 
        '\n'; 
    return false;
  }
  // Scan lines of the timestamps file for the EOH sentinel
  timestamps_ifstream.seekg(0);
  std::string line;
  std::string eoh_sentinel("%%ENDHEADER");
  std::streampos eoh_pos = -1;
  while (!timestamps_ifstream.eof()) {
    std::getline(timestamps_ifstream, line, '\n');
    // '%%ENDHEADER\n' must appear as the last characters of the header
    if ((line.find(eoh_sentinel) != std::string::npos) &&
        (eoh_sentinel.length() == line.length())) {
      eoh_pos = timestamps_ifstream.tellg();
      break;
    }
  }
  if (eoh_pos < 0) {
    std::cerr << "Timestamps file " << timestamps_filename << 
        " is missing end-of-header line " << eoh_sentinel << '\n';
    return false;
  }
  // Re-open timetamps file as binary to read the data
  timestamps_ifstream.close();
  timestamps_ifstream.open(timestamps_filename.c_str(), std::ifstream::binary);
  // Estimate the size of the data block
  timestamps_ifstream.seekg(0, std::ifstream::end);
  std::streampos data_size = timestamps_ifstream.tellg() - eoh_pos;
  // Seek to the putative start of the binary data block after header
  timestamps_ifstream.seekg(eoh_pos);
  if (timestamps_ifstream.rdstate()) {
    std::cerr << "Can not seek to position " << eoh_pos << " in " << 
        timestamps_filename << '\n';
    return false;
  }
  // Resize vector to match the number of (complete) timestamps
  timestamps_.resize(data_size/sizeof(uint32_t));
  // Read bytes from file and put them at &timestamps_ (is this dangerous?)
  if ((sizeof(uint32_t)*timestamps_.size() == data_size) && 
      (data_size > 0)) {
    // (is there a more elegant way to do this with iterator?)
    timestamps_ifstream.read((char*) &timestamps_[0], data_size); 
  }
  timestamps_ifstream.close();
  // Check that timestamps are monotonically increasing
  std::vector<uint32_t>::iterator timestamps_iter = timestamps_.begin();
  while (timestamps_iter != timestamps_.end()) {
    if (*(timestamps_iter+1) - *(timestamps_iter) <= 0) {
      std::cerr << "Timestamps violate monotonic increasing order\n";
      return false;
    }
    timestamps_iter++;
  }
  /* std::cout << "initialized timestamps\n"; */
  return true;
}


const bool& VideoReader::ready() const {
  return ready_;
}


const unsigned int& VideoReader::width() const {
  return width_;
}


const unsigned int& VideoReader::height() const {
  return height_;
}


const unsigned int& VideoReader::aspect_ratio_numerator() const {
  return aspect_ratio_numerator_;
}


const unsigned int& VideoReader::aspect_ratio_denominator() const {
  return aspect_ratio_denominator_;
}


const mpeg2convert_rgb_order_t& VideoReader::rgb_order() const {
  return RGB_ORDER;
}


const unsigned int& VideoReader::bits_per_pixel() const {
  return BITS_PER_PIXEL;
}


const framecount_t& VideoReader::number_of_frames() const {
  return number_of_frames_;
}


const framecount_t& VideoReader::current_index() const {
  return current_index_;
}


const uint32_t& VideoReader::current_timestamp() const {
  return current_timestamp_;
}


const std::vector<uint32_t>& VideoReader::all_timestamps() const {
  return timestamps_;
}


const uint8_t* VideoReader::current_picture() const {
  return current_picture_;
}


bool VideoReader::seek_frame(frameoff_t frame_off, std::ios::seekdir dir) {
  /* Check that mpeg_stream_ is open (it ought to remain persistently open until
  the destructor is called */
  if (!mpeg_stream_.is_open()) {
    std::cerr << "mpeg stream is not readable\n";
    return false;
  }
  /* It is important to keep the old values of current_picture_, current_index_,
  and current_timestamp_, so that the object reverts to original state in case
  new_index turns out to be invalid */
  framecount_t new_index;
  switch (dir) {
    case std::ios::beg: {
      new_index = 0 + frame_off;
      break;
    }
    case std::ios::cur: {
      new_index = current_index_ + frame_off;
      break;
    }
    case std::ios::end: {
      new_index = number_of_frames_ + frame_off;
      break;
    }
    default: {
      std::cerr << "Invalid seekdir flag\n";
      return false;
    }
  }
  if ((new_index >= 0) && (new_index < number_of_frames_)) {
    /* It is very important to call the clear() method, because otherwise
    lingerwise bad bits can impede the seekg! */
    mpeg_stream_.clear(); 
    mpeg_stream_.seekg(mpeg_index_[new_index]);
    /*
    std::cout << "\nnow at position " << mpeg_stream_.tellg() <<
        " in file, intended position is " << mpeg_index_[new_index] << 
        "\n";
    */
  } else {
    std::cerr << "Frame number " << new_index << " is out of bounds\n";
    return false;
  }

  /* Reset mpeg_decoder_ (we set the full_reset argument to 0 because we only
  seek within the same mpeg sequence) */
  mpeg2_reset(mpeg_decoder_, 0);
  /* std::cout << "reset mpeg decoder\n"; */

  uint8_t buffer[MPEG_BUFFER_SIZE];
  mpeg2_state_t state;
  bool completed_frame = 0;
  while (!completed_frame) {
    if (mpeg_stream_.rdstate() == std::ios::eofbit) {
      std::cerr << "reached end of mpeg file\n";
      return false;
    }
    if ((mpeg_stream_.rdstate() == std::ios::badbit) ||
        (mpeg_stream_.rdstate() == std::ios::failbit)) {
      std::cerr << "error reading mpeg file\n";
      return false;
    }
    state = mpeg2_parse(mpeg_decoder_);
    /* A quick tutorial on the simplified MPEG file format used by NSpike:

    A picture is a single captured video frame. All pictures in the mpeg files
    encoded with NSpike are I-frames, which means that they can be decoded as
    stand-alone images without reference to other video frames. A picture is
    further decomposed into one or more slices, which can be decoded
    individually and then assembled into the complete picture.

    Two or more pictures in sequence are packaged into a group-of-pictures
    (GOP).

    A video sequence is an ordered set of GOPs. A sequence starts with a
    sequence header and ends with an end-of-sequence code. The mpeg files that
    are saved by NSpike have a single sequence header at the beginning, which is
    the entry point for subsequence decoding.
    */
    switch (state) {
      case STATE_BUFFER: {
        // buffer is empty
        /* std::cout << "parser is ready for more data\n"; */
        mpeg_stream_.read(reinterpret_cast<char*>(buffer), sizeof(buffer));
        mpeg2_buffer(mpeg_decoder_, buffer, buffer + MPEG_BUFFER_SIZE);
        break;
      }
      case STATE_SEQUENCE: {
        /* std::cout << "sequence header was found\n"; */
        // set up output buffers to be in the desired pixel/color format
        mpeg2_convert(mpeg_decoder_, PIXEL_CONVERSION_SYMBOL, 0);
        break;
      }
      case STATE_SEQUENCE_REPEATED: {
        /* std::cout << "a repeated sequence header was found\n"; */
      }
      case STATE_GOP: {
        /* std::cout << "GOP header was found\n"; */
        break;
      }
      case STATE_PICTURE: {
        /* std::cout << "picture header was found\n"; */
        break;
      }
      case STATE_SLICE_1ST: {
        /* std::cout << "first slice of a multi-field picture was decoded\n"; */
        break;
      }
      case STATE_PICTURE_2ND: {
        /* std::cout << "a second field picture header was found\n"; */
        break;
      }
      case STATE_SLICE: {
        /* std::cout << "final slice of picture was decoded\n"; */
        break;
      }
      case STATE_END: {
        /* std::cout << "reached end of stream\n"; */
        return false;
        break;
      }
      case STATE_INVALID: {
        std::cerr << "buffer was filled without finding any codes\n";
        return false;
        break;
      }
      case STATE_INVALID_END: {
        std::cerr << "reached end of stream; unexpected termination\n";
        return false;
        break;
      }
      default: {
        std::cerr << "mpeg decoder buffer has unexpected state: " << 
            " mpeg_state_t enum value " << state << 
            " (see mpeg2.h to interpret)\n";
        return false;
      }
    }
    /* Did the decoder successfully populate the frame buffer? */
    if ((mpeg_info_->display_fbuf) && (mpeg_info_->display_fbuf->buf)) {
      /* display_fbuf->buf is a 3-element array of pointers to uint8_t arrays.
      When the frame is decoded as packed 8-bit RGB, only the first of these
      arrays is populated; display_fbuf->buf[1] and display_fbuf[2] are null
      pointers (reserved for YCbCr read-out?). display_fbuf->buf[0] contains
      width*height*3 uint8_t elements, ordered as follows: 

        R[row0,col0], G[row0,col0], B[row0,col0],
        R[row0,col1], G[row0,col1], B[row0,col1],
        R[row0,col2], G[row0,col2], B[row0,col2], 
        ...
        R[row1,col0], G[row1,col0], B[row1,col0],
        ...
        etc.  
      */
      current_picture_ = mpeg_info_->display_fbuf->buf[0];
      /* now we can safely replace the old value of current_index_ with
      new_index */
      current_index_ = new_index;
      current_timestamp_ = timestamps_[current_index_];
      completed_frame = true;
      /*
      std::cout << "frame number " << current_index_ <<
          " is ready (mpeg2_state_t enum value is " << state << ")\n";
      */
      return true;
    } 
  }
}


bool VideoReader::seek_frame(uint32_t target_timestamp) { 
  /* This is an overloaded version, which uses binary search to determine a
  frame_off argument to feed to the basic version of the method above */
  if (target_timestamp < timestamps_.front()) {
    std::cerr << "specified timestamp " << target_timestamp << 
        " is earlier than the first available timestamp (" << 
        timestamps_.front() << ") in the data\n";
    return false;
  }
  if (target_timestamp > timestamps_.back()) {
    std::cerr << "specified timestamp " << target_timestamp << 
        " is later than the last available timestamp (" << 
        timestamps_.back() << ") in the data\n";
    return false;
  }
  framecount_t new_index;
  frameoff_t low_idx = 0;
  frameoff_t high_idx = number_of_frames_ - 1;
  while (high_idx - low_idx > 1) {
    uint32_t low_ts = timestamps_[low_idx];
    uint32_t high_ts = timestamps_[high_idx];
    // Linear interpolation (remember that casting as unsigned integral type
    // truncates the fractional part of a float)
    new_index = static_cast<framecount_t> (static_cast<float> (low_idx) + 
        static_cast<float> (high_idx - low_idx) * 
        static_cast<float> (target_timestamp - low_ts) / 
        static_cast<float> (high_ts - low_ts));
    if (timestamps_[new_index] > target_timestamp) {
      high_idx = new_index;
    } else if (timestamps_[new_index] < target_timestamp) {
      low_idx = new_index;
    } else { // we found an exact match!
      high_idx = new_index;
      low_idx = new_index;
    }
  }
  // Call seek_frame again with an explicit frame_off argument
  if ((timestamps_[high_idx] - target_timestamp) < 
        (target_timestamp - timestamps_[low_idx])) {
    return seek_frame(high_idx, std::ios::beg);
  } else {
    // Note that this also covers the case of low_idx == high_idx
    return seek_frame(low_idx, std::ios::beg);
  }
}



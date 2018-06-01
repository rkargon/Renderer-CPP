//
//  common.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "common.h"

std::string containing_directory(const std::string &filename,
                                 const char *separator) {
  auto pos = filename.rfind(separator);
  if (pos == std::string::npos) {
    return ".";
  } else {
    return filename.substr(0, pos + 1);
  }
}

#ifndef __StringSplit_h_wanghan__
#define __StringSplit_h_wanghan__

#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>

namespace StringOperation{
  void split (const std::string & in,
	      std::vector<std::string > & out);
}

#endif

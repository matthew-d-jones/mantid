#ifndef MANTID_ISISREFLECTOMETRY_SEARCHRESULT_H
#define MANTID_ISISREFLECTOMETRY_SEARCHRESULT_H
#include <string>

namespace MantidQt {
namespace CustomInterfaces {

struct SearchResult {
  SearchResult() {}
  SearchResult(const std::string &runNumber, const std::string &desc,
               const std::string &loc)
      : runNumber(runNumber), description(desc), location(loc) {}
  std::string runNumber;
  std::string description;
  std::string location;
  std::string issues;
};
}
}
#endif // MANTID_ISISREFLECTOMETRY_SEARCHRESULT_H

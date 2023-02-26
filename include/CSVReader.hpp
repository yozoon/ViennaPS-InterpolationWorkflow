#ifndef PS_CSV_READER_HPP
#define PS_CSV_READER_HPP

#include <array>
#include <fstream>
#include <optional>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <lsSmartPointer.hpp>

#include "Utils.hpp"

template <class NumericType> class CSVReader {
  using ItemType = std::vector<NumericType>;
  using VectorType = std::vector<ItemType>;
  using ConstPtr = lsSmartPointer<const VectorType>;

  // Regex to find trailing and leading whitespaces
  const std::regex wsRegex = std::regex("^ +| +$|( ) +");

  std::string filename;
  char delimiter = ',';
  int numCols = 0;

  // An in-memory copy of the data
  VectorType data;
  std::unordered_map<std::string, NumericType> namedParameters;
  std::vector<NumericType> positionalParameters;

public:
  CSVReader() {}
  CSVReader(std::string passedFilename) : filename(passedFilename) {}

  void setFilename(std::string passedFilename) { filename = passedFilename; }

  // Returns a smart pointer to the in-memory copy of the data. If the in-memory
  // copy of the data is empty, the read function is called on the data source.
  ConstPtr getData() {
    if (data.empty())
      data = read();

    return ConstPtr::New(data);
  };

  std::vector<NumericType> getPositionalParameters() const {
    return positionalParameters;
  }

  std::unordered_map<std::string, NumericType> getNamedParameters() const {
    return namedParameters;
  }

private:
  void
  processPositionalParam(const std::string &input,
                         std::vector<NumericType> &positionalParameters) const {
    // Positional parameter
    auto v = psUtils::safeConvert<NumericType>(input);
    if (v.has_value())
      positionalParameters.push_back(v.value());
    else {
      std::cout << "Error while converting parameter '" << input
                << "' to numeric type.\n";
    }
  }

  void processNamedParam(
      const std::string &input,
      std::unordered_map<std::string, NumericType> &namedParameters) const {
    const std::string keyValueRegex =
        R"rgx(^[ \t]*([0-9a-zA-Z_-]+)[ \t]*=[ \t]*([0-9a-zA-Z_\-\.+]+)[ \t]*$)rgx";
    const std::regex rgx(keyValueRegex);

    std::smatch smatch;
    if (std::regex_search(input, smatch, rgx) && smatch.size() == 3) {
      auto v = psUtils::safeConvert<NumericType>(smatch[2]);
      if (v.has_value())
        namedParameters.insert({smatch[1], v.value()});
      else {
        std::cout << "Error while converting value of parameter '" << smatch[1]
                  << "'\n";
      }
    } else {
      std::cout << "Error while parsing parameter line '" << input << "'\n";
    }
  }

  void processParamLine(
      const std::string &line, std::vector<NumericType> &positionalParameters,
      std::unordered_map<std::string, NumericType> &namedParameters) const {
    std::istringstream iss(line);
    std::string tmp;

    // Skip the '#!' characters
    char c;
    iss >> c >> c;

    // Split the string at commas
    while (std::getline(iss, tmp, ',')) {
      if (tmp.find('=') == std::string::npos) {
        processPositionalParam(tmp, positionalParameters);
      } else {
        processNamedParam(tmp, namedParameters);
      }
    }
  }

  void processHeader() {
    auto opt = readHeader();
    if (opt.has_value()) {
      std::istringstream hdr(opt.value());
      std::string line;

      // Go over each comment line
      while (std::getline(hdr, line)) {
        // Check if the line is marked as a parameter line
        if (line.rfind("#!") == 0) {
          processParamLine(line, positionalParameters, namedParameters);
        }
      }
    }
  }

  std::optional<std::string> readHeader() {
    std::ifstream file(filename);
    std::string header;
    if (file.is_open()) {
      std::string line;
      // Iterate over each line
      while (std::getline(file, line)) {
        // Remove trailing and leading whitespaces
        line = std::regex_replace(line, wsRegex, "$1");

        // Skip empty lines at the top of the file
        if (line.empty())
          continue;

        // If the line is marked as comment and it is located before any data
        // at the top of the file, add it to the header string. Otherwise return
        // the header string, since we are now reading data.
        if (line.rfind('#') == 0) {
          header += '\n' + line;
        } else {
          file.close();
          break;
        }
      }
    } else {
      std::cout << "Couldn't open file '" << filename << "'\n";
      return {};
    }
    return {header};
  }

  std::optional<VectorType> readContent() {
    std::ifstream file(filename);
    if (file.is_open()) {
      auto data = VectorType{};

      std::string line;
      int lineCount = 0;

      // Iterate over each line
      while (std::getline(file, line)) {
        ++lineCount;

        // Remove trailing and leading whitespaces
        line = std::regex_replace(line, wsRegex, "$1");
        // Skip this line if it is marked as a comment
        if (line.rfind('#') == 0)
          continue;

        std::istringstream iss(line);
        std::string tmp;
        std::vector<NumericType> a;
        int i = 0;
        while (std::getline(iss, tmp, delimiter)) {
          auto valueOpt = psUtils::safeConvert<NumericType>(tmp);
          if (valueOpt)
            a.push_back(valueOpt.value());
          else {
            std::cout << "Error while reading line " << lineCount - 1 << " in '"
                      << filename << "'\n";
            return {};
          }
          ++i;
        }

        // The first row of actual data determins the data dimension
        if (numCols == 0)
          numCols = i;

        if (i != numCols) {
          std::cout << "Invalid number of columns in line " << lineCount - 1
                    << " in '" << filename << "'\n";
          return {};
        }
        data.push_back(a);
      }
      file.close();
      return data;
    } else {
      std::cout << "Couldn't open file '" << filename << "'\n";
      return {};
    }
  }

  VectorType read() {
    processHeader();
    auto contentOpt = readContent();
    if (contentOpt)
      return contentOpt.value();
    return {};
  }
};

#endif
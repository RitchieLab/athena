#ifndef STRINGMANIP_H_
#define STRINGMANIP_H_

#include <string>
#include <vector>

namespace data_manage
{

class Stringmanip
{
public:
    Stringmanip();
    /// check if number contained in string
    static int is_number(std::string str);
    /// convert to uppercase
    static std::string to_upper(std::string convert);
    /// splits string
    static int split_string(const std::string& input,
      std::string delim1, std::string delim2, std::vector<std::string>& results);
    static int stoi(std::string number);
    static unsigned stouint(std::string number);
    static int stodata(std::string number);
    static std::string itos(long number);
    static std::string itos(float number);
    static std::string itos(int number);
    static std::string itos(unsigned int number);
    static std::string itos(double number);
    static double stodouble(std::string number);
    static bool check_true_false(std::string value);
};

}

#endif /*STRINGMANIP_H_*/

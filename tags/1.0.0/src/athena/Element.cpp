//Element.cpp

#include "Element.h"

///
/// Constructor for types that contain a value
/// @param symbol string containing symbol used for the element
/// @param terminal_status bool that indicates whether this is a terminal element
/// @param priority_level int specifying how to prioritize evaluation for this operator
///
Element::Element(string symbol, bool terminal_status, int priority_level){
  name = symbol;
  isterminal = terminal_status;
  var_args = false;
  priority = priority_level;
  label = "nonterminal";
  shape = "box";
  style = "bold";
  type = label;
}

///
/// Constructor for types that contain a value
/// @param value 
/// @param terminal_status bool that indicates whether this is a terminal element
/// @param priority_level int specifying how to priority evaluation for this operator
///
Element::Element(DATATYPE value, bool terminal_status, int priority_level){
  stringstream ss;
  ss << value;
  name = ss.str();
  isterminal = terminal_status;
  var_args = false;
  priority = priority_level;
  label = "nonterminal";
  shape = "box";
  style = "bold";
  type = label;
}


///
/// Tracks whether Element is a terminal
/// If not a terminal then element must expand in the
/// grammar
/// @return terminal status
///
bool Element::get_isterminal() const{
  return isterminal;
}

///
/// Returns number of elements needed
/// for this element
/// @return number of elements
///
int Element::get_num_args() const{
  return num_args;
}

///
/// Sets the number of arguments
/// @param nargs number of arguments needed by the argument
///
void Element::set_num_args(int nargs){
  num_args = nargs;
}

///
/// Returns the symbol of the element
/// @return symbol 
///
string Element::get_name() const{
  return name;
}

///
/// Returns priority level for the element
/// @returns priority
///
int Element::get_priority() const{
  return priority;
}


///
/// Results are scaled from -1 to 1
/// @param x DATATYPE with result
/// @return scaled result
///
DATATYPE Element::ActivateSigmoid(DATATYPE x)
{
        if(x < -709) return -1.0;
        if(x > 709)  return 1.0;
        return(1.0 / (1.0 + exp(-x)));
}

///
/// Adjusts the result if it is infinite 
/// or not a number
/// @param x DATATYPE representing result
/// @return adjusted result 
///
DATATYPE Element::AdjustResult(DATATYPE x)
{
        if(isinf(x) == 1) {
                return 1.0;
        }
        if(isinf(x) == -1) {
                return -1.0;
        }
        if(isnan(x)) {
                return 0.0;
        }

        return x;
}

///
/// Overloads output for element
/// @param os ostream to write to
/// @param el Element to output
/// @returns output stream
///
ostream & operator << (ostream & os, const Element & el){
        os << el.name;
        return os;
}


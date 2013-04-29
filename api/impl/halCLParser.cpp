/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <sstream>
#include <iostream>
#include "halCLParser.h"
#include "hal.h"

using namespace std;
using namespace hal;

size_t CLParser::lineWidth = 85;

CLParser::CLParser() :
  _prefix("--"),
  _maxArgLen(0),
  _maxOptLen(0)
{
  addOptionFlag("help", "dsiplay this help page", false);
}

CLParser::~CLParser()
{
}
   
void CLParser::setOptionPrefix(const string& prefix)
{
  _prefix = prefix;
}

bool CLParser::hasOption(const string& name) const
{
  map<string, Option>::const_iterator i = _options.find(name);
  return i != _options.end() && i->second._flag == false;
}
   
void CLParser::addArgument(const string& name,
                           const string& description)
{
  if (hasOption(name) || hasFlag(name) || hasArgument(name))
  {
    throw std::runtime_error(string("name ") + name + "already present");
  }
  Argument arg = {name, description, ""};
  _args.push_back(arg);
  _maxArgLen = max(_maxArgLen, name.length());
}

bool CLParser::hasArgument(const string& name) const
{
  for (size_t i = 0; i < _args.size(); ++i)
  {
    if (_args[i]._name == name)
    {
      return true;
    }
  }
  return false;
}
   
void CLParser::addOptionFlag(const string& name,
                             const string& description,
                             bool defaultValue)
{
  if (hasOption(name) || hasFlag(name) || hasArgument(name))
  {
    throw std::runtime_error(string("name ") + name + "already present");
  }
  stringstream ss;
  ss << defaultValue;
  Option opt = { description, ss.str(), ss.str(), true };
  map<string, Option>::iterator i = _options.find(name);
  _options.insert(pair<string, Option>(name, opt));
  _maxOptLen = max(_maxOptLen, name.length());
}
   
bool CLParser::getFlag(const string& name) const
{
  bool value;
  map<string, Option>::const_iterator i = _options.find(name);
  if (i == _options.end() || i->second._flag == false)
  {
    throw hal_exception(string("Option flag") + name + " not recognized");
  }
  stringstream ss;
  try
  {
    ss << i->second._value;
    ss >> value;
  }
  catch (...)
  {
    throw hal_exception(string("type conversion getting flag ") + name);
  }
  return value;
}
   
bool CLParser::hasFlag(const string& name) const
{
  map<string, Option>::const_iterator i = _options.find(name);
  return i != _options.end() && i->second._flag == true;
}
   
void CLParser::parseOptions(int argc, char** argv)
{
  if (argc == 0)
  {
    return;
  }
  _exeName = argv[0];
  size_t argNum = 0;
  for (int i = 1; i < argc; ++i)
  {
    string name = argv[i];
    if (name.find(_prefix) == 0)
    {
      name = name.substr(_prefix.length());
      if (name == "help")
      {
        throw hal_exception(" ");
      }
      if (hasOption(name))
      {
        if (i == argc - 1)
        {
          throw hal_exception(string("Missing value for option ") + name);
        }
        map<string, Option>::iterator it = _options.find(name);
        it->second._value = argv[++i];
      }
      else if (hasFlag(name))
      {
        map<string, Option>::iterator it = _options.find(name);
        stringstream ss;
        ss << it->second._defaultValue;
        bool oldValue;
        ss >> oldValue;
        stringstream ss2;
        ss2 << !oldValue;
        it->second._value = ss2.str();
      }
      else
      {
        throw hal_exception(string("Unrecognized option: ") + string(argv[i]));
      }
    }
    else
    {
      if (argNum >= _args.size())
      {
        stringstream ss;
        ss << "Too many (" << argNum + 1 << ") arguments reached: " << name;
        throw hal_exception(ss.str());
      }
      _args[argNum++]._value = name;
    }
  }
  if (argNum != _args.size())
  {
    throw hal_exception("Too few (required positional) arguments");
  }
}

void CLParser::printUsage(ostream& os) const
{
  os << endl;
  if (_description.empty() == false)
  {
    stringstream vstring;
    vstring << HAL_VERSION;
    os << _exeName << " v" << vstring.str() << ": " 
       << multiLine(_description, _exeName.length() + 
                    vstring.str().length() + 4)
       << "\n\n";
  }
  os << "USAGE:\n" << _exeName << " [Options]";
  for (size_t i = 0; i < _args.size(); ++i)
  {
    os << " <" << _args[i]._name << ">";
  }
  os << endl << endl;
  os << "ARGUMENTS:\n";
  for (size_t i = 0; i < _args.size(); ++i)
  {
    string spacer(_maxArgLen - _args[i]._name.length(), ' ');
    os << _args[i]._name << ":   " << spacer
       << multiLine(_args[i]._description, _maxArgLen + 4) 
       << endl;
  }
  os << endl;
  os << "OPTIONS:\n";
  for (map<string, Option>::const_iterator i = _options.begin(); 
       i != _options.end(); ++i)
  {
    os << _prefix << i->first;
    size_t spacerLen = 8 + _maxOptLen - i->first.length();
    if (i->second._flag == false)
    {
      os << " <value>";
      spacerLen -= 8;
    }
    string spacer(spacerLen, ' ');
    os << ":   " << spacer  
       << multiLine(i->second._description + " [default = " +
                    i->second._defaultValue + "]", 14 + _maxOptLen) 
       << endl;
  }
  os << endl;
  if (_example.empty() == false)
  {
    os << "EXAMPLE:\n" << _exeName << " _example" << endl;
  }
  os << endl;
}

void CLParser::setDescription(const std::string& description)
{
  _description = description;
}

void CLParser::setExample(const std::string& example)
{
  _example = example;
}
   
string CLParser::multiLine(const string& line, size_t indent)
{
  if (indent >= lineWidth || line.length() + indent <= lineWidth)
  {
    return line;
  }
  string output;
  size_t width = 0;
  for (size_t i = 0; i < line.length(); ++i)
  {
    output += line[i];
    ++width;
    if (i != 0 && i != line.length() - 1 && width >= (lineWidth - indent))
    {
      bool inserted = false;
      if (!isspace(line[i]) && !isspace(line[i+1]))
      {
        for (size_t j = 1; j < output.length() && !inserted; ++j)
        {
          assert(i > j);
          if (isspace(line[i - j]))
          {
            inserted = true;
            output.insert((output.length() - j), 
                          string("\n") + string(indent, ' '));
            width = j;
          }
        }
      }
      if (!inserted)
      {
        output += '\n';
        output += string(indent, ' ');
        width = 0;
      }
    }    
  }
  return output;
}

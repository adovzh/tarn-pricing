#include <iostream>
#include <fstream>
#include <strstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include "Config.h"
#include "Logging.h"

namespace {
	// local handy type definitions
	typedef std::vector<std::string> StrVec;
	typedef std::map<std::string, StrVec> StrVecMap;
	typedef boost::shared_ptr<StrVecMap> StrVecMapPtr;

	// local auxiliary functions

	inline void push_element(StrVec& line, std::string& current)
	{
		boost::algorithm::trim(current);
		line.push_back(current);
		current.clear();
	}

	inline void push_line(StrVecMapPtr& dataMap, StrVec& line, std::string& current)
	{
		if (line.size() > 0)
		{
			push_element(line, current);
			std::string key;
			std::vector<std::string> values;
			bool first = true;

			for (std::vector<std::string>::const_iterator it = line.cbegin(); it != line.cend(); ++it)
			{
				if (first) { key = *it; first = false; }
				else values.push_back(*it);
			}

			dataMap->insert(std::pair<std::string, std::vector<std::string> >(key, values));
		}

		line.clear();
		current.clear();
	}
}

namespace tarnpricing {

int Config::getInt(const std::string& key) const
{
	StrVecMap::const_iterator it = data->find(key);

	if (it != data->cend()) return boost::lexical_cast<int>(it->second.front().c_str());
	else throw std::logic_error("No key found: " + key);
}

double Config::getDouble(const std::string& key) const
{
	StrVecMap::const_iterator it = data->find(key);

	if (it != data->cend()) return boost::lexical_cast<double>(it->second.front().c_str());
	else throw std::logic_error("No key found: " + key);
}

std::string Config::getString(const std::string& key) const
{
	StrVecMap::const_iterator it = data->find(key);

	if (it != data->cend()) return it->second.front();
	else throw std::logic_error("No key found: " + key);
}

void Config::getDoubles(const std::string& key, std::vector<double>& values) const
{
	StrVecMap::iterator it = data->find(key);

	if (it != data->end())
	{
		StrVec& strings = it->second;

		for (StrVec::const_iterator sit = strings.cbegin(); sit != strings.cend(); ++sit)
			values.push_back(boost::lexical_cast<double>(sit->c_str()));
	}
	else throw std::logic_error("No key found: " + key);
}

Config::ConstPtr Config::parse(const char* fileName)
{
	Config::ConstPtr ptr;
	std::ifstream f;	

	try
	{
		std::ios_base::iostate mask = f.exceptions();
		f.exceptions(mask | std::ios::failbit);
		f.open(fileName);
		f.exceptions(mask);
	}
	catch(std::ios_base::failure&)
	{
#ifdef ERROR_ENABLED
		const int bufferSize = 256;
		char buffer[bufferSize];
		std::strstream msg;
		msg << "Config::parse(\"" << fileName << "\")";
		if (!_strerror_s((char*)&buffer, bufferSize, msg.str())) ERROR_MESSAGE(buffer)
#endif

		return ptr;
	}

	int state = 0;
	int c = f.get();
	StrVec line;
	std::string current;
	StrVecMapPtr data(new StrVecMap);

	while (!f.eof())
	{
		switch(state)
		{
		case 0:
			switch(c)
			{
			case ',':
				push_element(line, current);
				break;
			case '\r': 
				break;
			case '\n':
				push_line(data, line, current);
				break;
			case '#':
				state = 1;
				break;
			default:
				current += c;
			}
			break;

		case 1:
			switch(c)
			{
			case '\n':
				f.putback(c);
				state = 0;
				break;
			}
		}

		c = f.get();
	}


	push_line(data, line, current);
	ptr.reset(new Config(data));

	return ptr;
}

} // namespace tarnpricing
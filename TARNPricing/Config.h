#ifndef __CONFIG_H
#define __CONFIG_H

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

namespace tarnpricing {

/**
 * Config class represents parsed information from a configuration CSV file. It can only be created 
 * via static parse() function, that loads it from a specified file. Once constructed, this class
 * provides operations to query config values by key. The key must exist and values must be of the correct type.
 *
 * std::logic_error is thrown if the key is absent. boost::bad_lexical_cast is thrown if the value is of the incorrect type.
 */
class Config
{
public:
	/// Retrieves the value from the configuration using the specified key and converts it to integer
	int getInt(const std::string& key) const;
	/// Retrieves the value from the configuration using the specified key and converts it to double
	double getDouble(const std::string& key) const;
	/// Retrieves the value from the configuration using the specified key
	std::string getString(const std::string& key) const;
	/// Retrieves the value from the configuration using the specified key and fills the specified vector with them.
	/// Values are converted to doubles
	void getDoubles(const std::string& key, std::vector<double>& values) const;

	typedef boost::shared_ptr<Config> Ptr;
	typedef boost::shared_ptr<const Config> ConstPtr;

	/// Parse the file, that fileName points to, create a fresh Config object and return a shared pointer to a constant object
	static ConstPtr parse(const char* fileName);
private:
	Config(const boost::shared_ptr<std::map<std::string, std::vector<std::string> > >& _data): data(_data) {}
	/// Internally config is just a map from keys (first column in CSV file) to a vector of strings
	/// Values are converted to a requested type on demand
	boost::shared_ptr<std::map<std::string, std::vector<std::string> > > data;
};

} // namespace tarnpricing

#endif // __CONFIG_H
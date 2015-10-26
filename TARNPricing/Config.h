#ifndef __CONFIG_H
#define __CONFIG_H

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

namespace tarnpricing {

class Config
{
public:
	int getInt(const std::string& key) const;
	double getDouble(const std::string& key) const;
	void getDoubles(const std::string& key, std::vector<double>& values) const;

	typedef boost::shared_ptr<Config> Ptr;
	typedef boost::shared_ptr<const Config> ConstPtr;

	static ConstPtr parse(const char* fileName);
private:
	Config(const boost::shared_ptr<std::map<std::string, std::vector<std::string> > >& _data): data(_data) {}
	boost::shared_ptr<std::map<std::string, std::vector<std::string> > > data;
};

} // namespace tarnpricing

#endif // __CONFIG_H
class AppCfg
{
public:
  AppCfg(std::ifstream &cfgfile, std::map<std::string, std::string> &options) 
  {
  }
  //AppCfg(const char* filename)
  //  : m_file(filename)
  //{
  //}
  AppCfg(const std::string& filename)
    : m_cfgfile(filename)
  {
    long begin, end;
    begin = m_cfgfile.tellg();
    m_cfgfile.seekg (0, std::ios::end);
    end = m_cfgfile.tellg();
    //m_cfgfile.close();
    m_cfgfile.seekg (0, std::ios::beg);
    std::cout << "DBG size: " << (end-begin) << " bytes." << std::endl;
    auto size = std::filesystem::file_size(filename);
    std::cout << "DBG size: " << size << " bytes." << std::endl;
	
  }

  //void parse(std::ifstream &cfgfile, std::map<std::string, std::string> &options)
  void parse(std::map<std::string, std::string> &options)
  {
    std::cout << "AppCfg parse " << std::endl;
    for (std::string line; getline(m_cfgfile, line); )
    {
        //std::cout << "DBG line " << line << std::endl;
	std::istringstream iss(line);
	//std::istringstream iss;
	//iss.str(line);
        //std::cout << "DBG line " << iss.str() << std::endl;
	std::string id, eq, val;

        bool error = false;

        if (!(iss >> id))
        {
            //std::cout << "DBG error: " << iss.str() << " to " << id << std::endl;
            //std::cout << "DBG error: " << iss.str() << std::endl;
            error = true;
        }
        else if (id[0] == '#')
        {
            //std::cout << "DBG id " << id[0] << " " << iss.str() << " ID " << id << std::endl;
            continue;
        }
        else if (!(iss >> eq >> val) || eq != "=" || iss.get() != EOF)
        {
            error = true;
        }

        if (error)
        {
            // do something appropriate: throw, skip, warn, etc.
	    //std::cout << "DBG error: " << line << std::endl;
        }
        else
        {
	    //std::cout << "DBG id val " << id << " " << val << std::endl;
            options[id] = val;
        }
    }
    m_cfg = options;
  }

  void print(std::map<std::string, std::string> &options)
  {
    std::cout << "AppCfg print " << std::endl;
    std::cout << " Options: " << std::endl;
    for (auto& t : options) {
      std::cout << t.first << " : " << t.second << std::endl;
    }
  }

private:
  std::ifstream m_cfgfile;
  std::map<std::string, std::string> m_cfg;
};


bool to_bool(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}


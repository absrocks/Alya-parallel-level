/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "string_dict.h"

#include "forwarding_logger.h"

#include <cassert>
#include <iterator>
#include <limits>

#include "tokenizer.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
// report error conversion in default.
bool String_dict::S_is_error_conversion_report_on = true;

//----------------------------------------------------------------------
String_dict::String_dict()
{
    // empty
}

//----------------------------------------------------------------------
String_dict::String_dict(String_dict const & dict)
{
    // iterate through all interfaces of dict
    // and insert them here.
    for ( const_iterator i = dict.m_map.begin(); i != dict.m_map.end(); ++i) {
        insert( i->first, i->second);
    }
}

//----------------------------------------------------------------------
String_dict::~String_dict()
{
    this->clear();
}

//----------------------------------------------------------------------
std::pair< String_dict::iterator, bool> String_dict::insert(
    std::string const & key,
    std::string const & value)
{
    std::pair< iterator, bool> ret = m_map.insert( value_type( key, value));
    if ( !(ret.second)) {
        // Key already in std::set does not change value, we do this explicitly.
        (ret.first)->second = value;
    }
    return ret;
}

//----------------------------------------------------------------------
mi::Sint32 String_dict::insert_all(String_dict const & dict)
{
    mi::Sint32 n_inserted = 0;
    for(const_iterator ii = dict.m_map.begin(); ii != dict.m_map.end(); ++ii) {
        std::pair< iterator, bool> p = this->insert(ii->first, ii->second);
        n_inserted += p.second ? 1 : 0;
    }
    return n_inserted;
}

//----------------------------------------------------------------------
mi::Sint32 String_dict::insert_new(String_dict const & dict)
{
    mi::Sint32 n_inserted = 0;
    for(const_iterator ii = dict.m_map.begin(); ii != dict.m_map.end(); ++ii) {
        if (m_map.find(ii->first) == m_map.end())
        {
            insert(ii->first, ii->second);
            n_inserted++;
        }
    }
    return n_inserted;
}

//----------------------------------------------------------------------
std::string String_dict::get(std::string const & key,
                             std::string const & default_value) const
{
    const_iterator pos = m_map.find( key);
    if (pos != m_map.end()) {
        return pos->second;
    }
    return default_value;
}

//----------------------------------------------------------------------
void String_dict::set(std::string const & key,
                            std::string const & value )
{
    this->insert(key, value);
}

//----------------------------------------------------------------------
void String_dict::clear()
{
    m_map.clear();
}

//----------------------------------------------------------------------
bool String_dict::is_defined(
    const std::string& key) const
{
    return (m_map.find(key) != m_map.end());
}

//----------------------------------------------------------------------
size_t String_dict::size() const
{
    return m_map.size();
}

//----------------------------------------------------------------------
String_dict::const_iterator String_dict::begin() const
{
    return m_map.begin();
}

//----------------------------------------------------------------------
String_dict::const_iterator String_dict::end() const
{
    return m_map.end();
}

//----------------------------------------------------------------------
void String_dict::erase(const std::string & key)
{
    iterator ii = m_map.find(key);
    if (ii != m_map.end()) {
        m_map.erase(ii);
    }
}

//----------------------------------------------------------------------
void String_dict::write(std::ostream & os,
                        std::string const & prefix,
                        bool is_sort) const
{
    if(is_sort){
        std::vector< std::string > sortbuf;
        for(const_iterator ii = this->begin(); ii != this->end(); ++ii){
            sortbuf.push_back(prefix + ii->first + " = " + ii->second);
        }
        std::sort(sortbuf.begin(), sortbuf.end());
        size_t const entry_count = sortbuf.size();
        for(size_t i = 0; i < entry_count; ++i){
            os << sortbuf[i] << "\n";
        }
    }
    else{
        for(const_iterator ii = this->begin(); ii != this->end(); ++ii){
            os << prefix << ii->first << " = " << ii->second << "\n";
        }
    }
}

//----------------------------------------------------------------------
void String_dict::read(std::istream & is)
{
    mi::Sint32 nline = 0;              // line number
    std::string line;
    while (getline(is, line, '\n'))
    {
        ++nline;
        if (line.empty())
        {
            continue;           // empty line
        }
        else if (line[0] == '#')
        {
            continue;           // comment line
        }
        else if (line.find_first_not_of(" \t") == std::string::npos)
        {
            continue;           // white space line
        }
        std::string::size_type pos_equal = line.find("=");
        std::string::size_type pos_multiline = line.find("<<");
        if (pos_equal == std::string::npos && pos_multiline == std::string::npos)
        {
            ERROR_LOG << "String_dict: Cannot find '=' or '<<', unrecognized entry "
                      << "ignored at line " << nline << ": [" << line << "].";
            continue; // unrecognized line
        }

        // Handle single line ("a=b") and multi-line entries ("a<<(END)\n...values...\n(END)")
        std::string::size_type pos = pos_equal;
        bool is_multline = false;
        if (pos_multiline != std::string::npos &&
            (pos_equal == std::string::npos || pos_multiline < pos_equal))
        {
            is_multline = true;
            pos = pos_multiline;
        }

        std::string key = line.substr(0, pos);
        std::string val = line.substr(pos + (is_multline ? 2 : 1));

        // Trim white space
        key.erase(0, key.find_first_not_of(" \n\r\t"));
        key.erase(key.find_last_not_of(" \n\r\t") + 1);

        val.erase(0, val.find_first_not_of(" \n\r\t"));
        val.erase(val.find_last_not_of(" \n\r\t") + 1);

        std::string multiline_delimiter;
        if (is_multline)
        {
            // Multi-line entries are similar to "here documents" in bash, meaning the delimiter
            // following the "<<" is specified by the user, and input will be read until a line
            // containing only the delimiter is found.
            if (val.empty())
            {
                ERROR_LOG << "Missing multi-line delimiter for key '"
                         << key << "' starting at line " << nline << ".";
                continue; // ignore this entry
            }

            // Take the value as the delimiter, the actual multi-line value is read below
            multiline_delimiter = val;
            val = "";

            // Read multi-line value, i.e. all lines until the delimiter is found
            bool delimiter_found = false;
            mi::Sint32 multiline_start = nline;
            while (!delimiter_found && getline(is, line, '\n'))
            {
                ++nline;

                // Find lines starting with the delimiter
                if (line.substr(0, multiline_delimiter.size()) == multiline_delimiter)
                {
                    // The must be nothing but  whitespace after the delimiter
                    if (line.find_first_not_of(" \n\r\t", multiline_delimiter.size()) == std::string::npos)
                    {
                        delimiter_found = true;
                        break;
                    }
                }

                val += line + "\n";
            }

            if (!delimiter_found)
            {
                ERROR_LOG << "Reached EOF while reading multi-line entry for key '" << key << "' "
                          << "starting at line " << multiline_start << ", delimiter is '" << multiline_delimiter << "'.";
                continue; // ignore this entry
            }
        }

        this->insert(key, val);
    }
}


//----------------------------------------------------------------------
bool String_dict::operator==(String_dict const & rhs ) const
{
    if ( rhs.m_map.size( ) != m_map.size( ) ){
        return false;
    }
    for ( const_iterator it = rhs.begin( ); it != rhs.end( ); ++it ) {
        const_iterator found = m_map.find( it->first );
        if ( found == m_map.end( ) ) {
            return false;
        }
        if ( found->second != it->second ) {
            return false;
        }
    }
    return true;
}

//----------------------------------------------------------------------
bool String_dict::operator!=(String_dict const & rhs ) const
{
    return !(*this==rhs);
}

//----------------------------------------------------------------------
bool is_all_keys_defined(String_dict const & option,
                         std::vector< std::string > const & key_list,
                         std::vector< std::string > * p_undef_list)
{
    size_t const list_count = key_list.size();
    bool all_defined = true;
    for(size_t i = 0; i < list_count; ++i){
        if(!option.is_defined(key_list[i])){
            all_defined = false;
            if(p_undef_list != 0){
                p_undef_list->push_back(key_list[i]);
            }
        }
    }

    return all_defined;
}

//----------------------------------------------------------------------
bool is_all_keys_defined(String_dict const & option,
                         char const * const p_key_list[],
                         std::vector< std::string > * p_undef_list)
{
    std::vector< std::string > key_list;
    for(mi::Sint32 i = 0; p_key_list[i] != 0; ++i){
        key_list.push_back(std::string(p_key_list[i]));
    }
    return is_all_keys_defined(option, key_list, p_undef_list);
}

//----------------------------------------------------------------------
bool is_all_keys_defined_with_string_report(
    String_dict const & option,
    std::vector< std::string > const & key_list,
    std::string & undef_list_str)
{
    std::vector< std::string > undef_list;
    if(!is_all_keys_defined(option, key_list, &undef_list)){
        std::stringstream sstr;
        std::copy(undef_list.begin(), undef_list.end(),
                  std::ostream_iterator< std::string >(sstr, " "));
        undef_list_str = sstr.str();
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
bool is_all_keys_defined_with_string_report(
    String_dict const & option,
    char const * const p_key_list[],
    std::string & undef_list_str)
{
    std::vector< std::string > undef_list;
    if(!is_all_keys_defined(option, p_key_list, &undef_list)){
        std::stringstream sstr;
        std::copy(undef_list.begin(), undef_list.end(),
                  std::ostream_iterator< std::string >(sstr, " "));
        undef_list_str = sstr.str();
        return false;
    }
    return true;
}


//----------------------------------------------------------------------
void string_array_to_string_dict(mi::Sint32  const    argc,
                                 char *        argv[],
                                 String_dict & option,
                                 bool set_true_if_key_only)
{
    if(argc >= 1){
        option.insert("command:", std::string(argv[0]));
        // keep the full command line also

        std::stringstream sstr;
        for(mi::Sint32 i = 0; i < argc; ++i){
            sstr << argv[i] << " ";
        }
        option.insert("command_line:", sstr.str());
    }

    mi::Sint32 arg_n = 0;
    mi::Sint32 i     = 1;
    while(i < argc){
        // case -key
        if(argv[i][0] == '-'){
            std::string key = "<empty_key>";
            if(argv[i][1] != '\0'){
                key = std::string(static_cast< const char * >(&(argv[i][1])));
            }
            ++i;
            std::string value = "<no_value>";
            if(set_true_if_key_only && option.is_defined(key)){
                value = option.get(key); // get default
            }
            if((i < argc) && (argv[i][0] != '-')){
                // a value found.
                value = argv[i];
                ++i;
            }else{
                if(set_true_if_key_only){
                    // value = "1" means bool true. When only option is
                    // defined, but no value found, this becomes true.
                    value = "1";
                }
            }
            option.insert(key, value);
        }
        else{
            std::stringstream sstr;
            sstr << "arg_" << arg_n << ":";
            option.insert(sstr.str(), std::string(argv[i]));
            ++arg_n;
            ++i;
        }
    }
}

//----------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os,
                          String_dict const & dict)
{
    for(String_dict::const_iterator ii = dict.begin(); ii != dict.end(); ++ii){
        os << ii->first << " = " << ii->second << "\n";
    }
    return os;
}


//----------------------------------------------------------------------
static bool get_file_magic_sub(std::string const & fname,
                               FILE * pfp,
                               std::string const & magic_string,
                               mi::Sint32 &file_version,
                               std::string & error_mes)
{
    assert(pfp != 0);

    const mi::Sint32 BUFSIZE = 1024;
    char line[BUFSIZE];
    if(fgets(line, BUFSIZE, pfp) == 0){
        error_mes = "Cannot read the first line of the file [" + fname + "]";
        return false;
    }

    std::string const first_line(line);
    std::string const should_magic_prefix = "#! " + magic_string;

    // first line should 'should_magic_prefix' + ' k' (space + one
    // digit) at least. The '+ 2' comes from this minimal two
    // characters ' k'.
    if(first_line.size() < (should_magic_prefix.size() + 2)){
        size_t len = first_line.find("\n"); // remove \n
        error_mes = std::string("Too short file magic [") + first_line.substr(0, len) + "], expected [" +
            should_magic_prefix + " {version_int}], " +
            "where {version_int} is a integer number of the project file version. ex. '0'.";
        return false;
    }

    std::string const magic_prefix = first_line.substr(0, should_magic_prefix.size());
    std::string const ver_str = first_line.substr(should_magic_prefix.size() + 1);

    //     std::cout << "should magic [" << should_magic_prefix << "]\n"
    //               << "see    magic [" << magic_prefix        << "]\n"
    //               << "ver str      [" << ver_str             << "]\n";

    if(magic_prefix != should_magic_prefix){
        error_mes = "Wrong file magic. It should start with [" +  should_magic_prefix + "]";
        return false;
    }

    mi::Sint32 ver_num = -1;
    mi::Sint32 const ret_elem = sscanf(ver_str.c_str(), "%6d", &ver_num);
    if(ret_elem != 1){
        error_mes = "Cannot find version number.";
        return false;
    }

    file_version = ver_num;

    return true;
}

//----------------------------------------------------------------------
mi::Sint32 get_file_magic(std::string const & fname,
                          std::string const & magic_string,
                          mi::Sint32 &        file_version,
                          std::string *       p_error_mes)
{
    FILE * pfp = fopen(fname.c_str(), "rb");
    if(pfp == 0){
        if(p_error_mes != 0){
            (*p_error_mes)= "Cannot open [" + fname + "]";
        }
        return -1;
    }

    std::string err_mes;
    bool const ret = get_file_magic_sub(fname, pfp, magic_string, file_version,
                                        err_mes);
    if(p_error_mes != 0){
        (*p_error_mes) = err_mes;
    }

    fclose(pfp);

    if(!ret){
        return -2;              // wrong header
    }

    return 0;                   // success
}

//----------------------------------------------------------------------
bool get_string_dict_with_magic(std::string const & dict_opt_fname,
                                std::string const & magic_str,
                                String_dict       & dict_opt,
                                mi::Sint32        & file_version,
                                std::string       * p_error_mes)
{
    mi::Sint32 fver = -1;
    const mi::Sint32 ret = get_file_magic(dict_opt_fname, magic_str, fver, p_error_mes);
    if((ret == -1) || (ret == -2)){
        // -1 ... file not found, -2 ... wrong header, error message is already in there
        return false;
    }
    assert(ret == 0);           // we know only two error types, others should be OK.

    // the magic is correct.
    file_version = fver;
    dict_opt.clear();

    // Replaced with stringstream since we can not use input file
    // stream by a compiler bug.
    std::string data_str, error_mes;
    bool const is_success = load_string_from_file(dict_opt_fname, data_str, error_mes);
    if(!is_success){
        ERROR_LOG << "get_string_dict_with_magic: " << error_mes;
    }
    assert(is_success);

    std::istringstream isstr(data_str);
    dict_opt.read(isstr);

    return true;
}

//----------------------------------------------------------------------
mi::Sint32 string_dict_key_prefix_filter(
    String_dict const & src_dict,
    std::string const & prefix,
    String_dict       & filtered,
    bool                is_delete_prefix)
{
    static std::string const fh = "dictionary_key_prefix_filter: ";
    if(prefix == ""){           // empty prefix
        ERROR_LOG << fh << "invalid prefix string (empty)";
        return -1;
    }

    mi::Sint32 n_found = 0;
    for(String_dict::const_iterator di = src_dict.begin(); di != src_dict.end(); ++di){
        if(di->first.find(prefix) == 0){
            std::string key = di->first;
            if(is_delete_prefix){
                if(key == prefix){
                    ERROR_LOG << fh << "key = prefix. Cannot remove the whole key";
                    return -1;
                }
                key = key.substr(prefix.size());
            }
            filtered.insert(key, di->second);

            ++n_found;
        }
    }

    return n_found;
}

//----------------------------------------------------------------------
// type conversion helpers

/// show conversion error report
///
/// \param[in] type_name conversion type name to (from string)
/// \param[in] value_str     input string that conversion failed
void string_conversion_error_report(std::string const & type_name, std::string const & value_str)
{
    if(String_dict::is_error_conversion_report_on()){
        ERROR_LOG << "String_dict conversion error: string to " << type_name
                  << " string = [" << value_str << "]";
    }
}

//----------------------------------------------------------------------
bool get_bool(std::string const & value_str)
{
    if (value_str == "1" || value_str == "true" || value_str == "yes" || value_str == "on") {
        return true;
    }
    if (value_str == "0" || value_str == "false" || value_str == "no" || value_str == "off") {
        return false;
    }
    string_conversion_error_report("bool", value_str);
    return false;
}

std::string to_string(bool value)
{
    if(value){
        return std::string("1");
    }
    return std::string("0");
}

//----------------------------------------------------------------------
mi::Sint32 get_sint32(std::string const & value_str)
{
    mi::Sint32 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Sint32", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Sint32 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}


//----------------------------------------------------------------------
mi::Uint32 get_uint32(std::string const & value_str)
{
    mi::Uint32 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Uint32", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Uint32 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

//----------------------------------------------------------------------
mi::Sint64 get_sint64(std::string const & value_str)
{
    mi::Sint64 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Sint64", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Sint64 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

//----------------------------------------------------------------------
mi::Uint64 get_uint64(std::string const & value_str)
{
    mi::Uint64 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Uint64", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Uint64 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

//----------------------------------------------------------------------
mi::Float32 get_float32(std::string const & value_str)
{
    mi::Float32 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Float32", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Float32 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}


//----------------------------------------------------------------------
mi::Float64 get_float64(std::string const & value_str)
{
    mi::Float64 ret = 0;
    std::istringstream sstr(value_str);
    sstr >> ret;
    if (sstr.fail()) {
        string_conversion_error_report("mi::Float64", value_str);
        ret = 0;
    }
    return ret;
}

std::string to_string(mi::Float64 value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Sint32, 2 > get_vec_sint32_2(std::string const & value_str)
{
    mi::Sint32 const DIM = 2;
    mi::math::Vector< mi::Sint32, DIM > ret(0,0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Sint32," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Sint32, DIM > zero = {0, 0,};
            return zero;
        }
    }
    mi::math::Vector_struct< mi::Sint32, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Sint32, 2 > const & value)
{
    mi::Sint32 const DIM = 2;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Sint64, 2 > get_vec_sint64_2(std::string const & value_str)
{
    mi::Sint64 const DIM = 2;
    mi::math::Vector< mi::Sint64, DIM > ret(0,0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Sint64," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Sint64, DIM > zero = {0, 0,};
            return zero;
        }
    }
    mi::math::Vector_struct< mi::Sint64, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Sint64, 2 > const & value)
{
    mi::Sint32 const DIM = 2;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Sint32, 3 > get_vec_sint32_3(std::string const & value_str)
{
    mi::Sint32 const DIM = 3;
    mi::math::Vector< mi::Sint32, DIM > ret(0, 0, 0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Sint32," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Sint32, DIM > zero = {0, 0, 0,};
            return zero;
        }
    }

    mi::math::Vector_struct< mi::Sint32, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    ret_st.z = ret.z;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Sint32, 3 > const & value)
{
    mi::Sint32 const DIM = 3;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}


//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Uint32, 2 > get_vec_uint32_2(std::string const & value_str)
{
    mi::Sint32 const DIM = 2;
    mi::math::Vector< mi::Uint32, DIM > ret(0,0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Uint32," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Uint32, DIM > zero = {0, 0,};
            return zero;
        }
    }
    mi::math::Vector_struct< mi::Uint32, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Uint32, 2 > const & value)
{
    mi::Sint32 const DIM = 2;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Uint32, 3 > get_vec_uint32_3(std::string const & value_str)
{
    mi::Sint32 const DIM = 3;
    mi::math::Vector< mi::Uint32, DIM > ret(0, 0, 0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Uint32," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Uint32, DIM > zero = {0, 0, 0,};
            return zero;
        }
    }

    mi::math::Vector_struct< mi::Uint32, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    ret_st.z = ret.z;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Uint32, 3 > const & value)
{
    mi::Sint32 const DIM = 3;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Sint64, 3 > get_vec_sint64_3(std::string const & value_str)
{
    mi::Sint32 const DIM = 3;
    mi::math::Vector< mi::Sint64, DIM > ret(0, 0, 0);
    std::istringstream isstr(value_str);
    for(mi::Sint32 i = 0; i < DIM; ++i){
        isstr >> ret[i];
        if(isstr.fail()) {
            std::stringstream sstr;
            sstr << "mi::math::Vector<Sint64," << DIM << ">";
            string_conversion_error_report(sstr.str(), value_str);
            mi::math::Vector_struct< mi::Sint64, DIM > zero = {0, 0, 0,};
            return zero;
        }
    }

    mi::math::Vector_struct< mi::Sint64, DIM > ret_st;
    ret_st.x = ret.x;
    ret_st.y = ret.y;
    ret_st.z = ret.z;
    return ret_st;
}

std::string to_string(mi::math::Vector< mi::Sint64, 3 > const & value)
{
    mi::Sint32 const DIM = 3;
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < DIM; ++i){
        sstr << value[i];
        if(i < (DIM - 1)){
            sstr << " ";
        }
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Float32, 2 > get_vec_float32_2(std::string const & value_str)
{
    mi::math::Vector< mi::Float32, 2 > ret(0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 2; ++i){
        sstr >> ret[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Vector<Float32,2>", value_str);
            return mi::math::Vector< mi::Float32, 2 >(0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Vector< mi::Float32, 2 > const & value)
{
    std::stringstream sstr;
    sstr << value[0] << " " << value[1];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Float32, 3 > get_vec_float32_3(std::string const & value_str)
{
    mi::math::Vector< mi::Float32, 3 > ret(0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Vector<Float32,3>", value_str);
            return mi::math::Vector< mi::Float32, 3 >(0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Vector< mi::Float32, 3 > const & value)
{
    std::stringstream sstr;
    sstr << value[0] << " " << value[1] << " " << value[2];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Float32, 4 > get_vec_float32_4(std::string const & value_str)
{
    mi::math::Vector< mi::Float32, 4 > ret(0.0f, 0.0f, 0.0f, 0.0f);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 4; ++i){
        sstr >> ret[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Vector<Float32,4>", value_str);
            return mi::math::Vector< mi::Float32, 4 >(0.0f, 0.0f, 0.0f, 0.0f);
        }
    }
    return ret;
}

std::string to_string(mi::math::Vector< mi::Float32, 4 > const & value)
{
    std::stringstream sstr;
    sstr << value[0] << " " << value[1] << " " << value[2] << " " << value[3];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Matrix_struct< mi::Float32, 4, 4 > get_mat_float32_4_4(
    std::string const & value_str)
{
    // Read the values in the same order as they are printed by the logger
    mi::math::Matrix< mi::Float32, 4, 4 > mat44;
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 4; ++i){ // n
        for(mi::Sint32 j = 0; j < 4; ++j){ // m
            mi::Float32 elem = 0.0;
            sstr >> elem;
            if(sstr.fail()) {
                string_conversion_error_report("mi::math::Matrix<Float32,4,4>", value_str);
                mi::math::Matrix< mi::Float32, 4, 4 > zero(0.f);
                return zero;
            }
            mat44.set(j, i, elem);
        }
    }

    return mat44;
}

std::string to_string(mi::math::Matrix< mi::Float32, 4, 4 > const & value)
{
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < 16; ++i){
        sstr << value.get(i) << " ";
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Matrix_struct< mi::Float64, 4, 4 > get_mat_float64_4_4(
    std::string const & value_str)
{
    mi::math::Matrix< mi::Float64, 4, 4 > mat44;
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 4; ++i){ // n
        for(mi::Sint32 j = 0; j < 4; ++j){ // m
            mi::Float64 elem = 0.0;
            sstr >> elem;
            if(sstr.fail()) {
                string_conversion_error_report("mi::math::Matrix<Float64,4,4>", value_str);
                mi::math::Matrix_struct< mi::Float64, 4, 4 > zero = {
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0,
                };
                return zero;
            }
            mat44.set(i, j, elem);
        }
    }

    // rely on operator=()
    mi::math::Matrix_struct< mi::Float64, 4, 4 > mat44_st = mat44;
    return mat44_st;
}

std::string to_string(mi::math::Matrix< mi::Float64, 4, 4 > const & value)
{
    std::stringstream sstr;
    for(mi::Sint32 i = 0; i < 16; ++i){
        sstr << value.get(i) << " ";
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Color_struct get_color(std::string const & value_str)
{
    mi::math::Color ret(0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 4; ++i){
        sstr >> ret[i];
        if(sstr.fail()) {
            // Replace missing alpha component with 1
            if (i == 3 && sstr.eof())
            {
                ret[i] = 1.f;
                break;
            }
            string_conversion_error_report("mi::math::Color", value_str);
            return mi::math::Color(0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Color const & value)
{
    std::stringstream sstr;
    sstr << value[0] << " " << value[1] << " " << value[2] << " " << value[3];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Float32, 2 > get_bbox_float32_2(std::string const & value_str)
{
    mi::math::Bbox< mi::Float32, 2 > ret(0.0f, 0.0f, 0.0f, 0.0f);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 2; ++i){
        sstr >> ret.min[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float32,2>", value_str);
            return mi::math::Bbox< mi::Float32, 2 >(0.0f, 0.0f, 0.0f, 0.0f);
        }
    }
    for(mi::Sint32 i = 0; i < 2; ++i){
        sstr >> ret.max[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float32,2>", value_str);
            return mi::math::Bbox< mi::Float32, 2 >(0.0f, 0.0f, 0.0f, 0.0f);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Float32, 2 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1]
         << value.max[0] << " " << value.max[1];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Float64, 2 > get_bbox_float64_2(std::string const & value_str)
{
    mi::math::Bbox< mi::Float64, 2 > ret(0.0, 0.0, 0.0, 0.0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 2; ++i){
        sstr >> ret.min[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float64,2>", value_str);
            return mi::math::Bbox< mi::Float64, 2 >(0.0, 0.0, 0.0, 0.0);
        }
    }
    for(mi::Sint32 i = 0; i < 2; ++i){
        sstr >> ret.max[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float64,2>", value_str);
            return mi::math::Bbox< mi::Float64, 2 >(0.0, 0.0, 0.0, 0.0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Float64, 2 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1]
         << value.max[0] << " " << value.max[1];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Float32, 3 > get_bbox_float32_3(std::string const & value_str)
{
    mi::math::Bbox< mi::Float32, 3 > ret(0, 0, 0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.min[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float32,3>", value_str);
            return mi::math::Bbox< mi::Float32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.max[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Float32,3>", value_str);
            return mi::math::Bbox< mi::Float32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Float32, 3 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1] << " " << value.min[2] << " "
         << value.max[0] << " " << value.max[1] << " " << value.max[2];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Uint32, 3 > get_bbox_uint32_3(std::string const & value_str)
{
    mi::math::Bbox< mi::Uint32, 3 > ret(0, 0, 0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.min[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Uint32,3>", value_str);
            return mi::math::Bbox< mi::Uint32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.max[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Uint32,3>", value_str);
            return mi::math::Bbox< mi::Uint32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Uint32, 3 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1] << " " << value.min[2] << " "
         << value.max[0] << " " << value.max[1] << " " << value.max[2];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Sint64, 3 > get_bbox_sint64_3(std::string const & value_str)
{
    mi::math::Bbox< mi::Sint64, 3 > ret(0, 0, 0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.min[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Sint64,3>", value_str);
            return mi::math::Bbox< mi::Sint64, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    for(mi::Sint32 i = 0; i < 3; ++i){
        sstr >> ret.max[i];
        if(sstr.fail()) {
            string_conversion_error_report("mi::math::Bbox<Sint64,3>", value_str);
            return mi::math::Bbox< mi::Sint64, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Sint64, 3 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1] << " " << value.min[2] << " "
         << value.max[0] << " " << value.max[1] << " " << value.max[2];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Sint32, 2 > get_bbox_sint32_2(std::string const & value_str)
{
    mi::math::Bbox< mi::Sint32, 2 > ret(0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 2; ++i)
    {
        sstr >> ret.min[i];
        if(sstr.fail())
        {
            string_conversion_error_report("mi::math::Bbox<Sint32,2>", value_str);
            return mi::math::Bbox< mi::Sint32, 2 >(0, 0, 0, 0);
        }
    }
    for(mi::Sint32 i = 0; i < 2; ++i)
    {
        sstr >> ret.max[i];
        if(sstr.fail())
        {
            string_conversion_error_report("mi::math::Bbox<Sint32,2>", value_str);
            return mi::math::Bbox< mi::Sint32, 2 >(0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Sint32, 2 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1] << " " 
         << value.max[0] << " " << value.max[1];
    return sstr.str();
}

//----------------------------------------------------------------------
mi::math::Bbox_struct< mi::Sint32, 3 > get_bbox_sint32_3(std::string const & value_str)
{
    mi::math::Bbox< mi::Sint32, 3 > ret(0, 0, 0, 0, 0, 0);
    std::istringstream sstr(value_str);
    for(mi::Sint32 i = 0; i < 3; ++i)
    {
        sstr >> ret.min[i];
        if(sstr.fail())
        {
            string_conversion_error_report("mi::math::Bbox<Sint32,3>", value_str);
            return mi::math::Bbox< mi::Sint32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    for(mi::Sint32 i = 0; i < 3; ++i)
    {
        sstr >> ret.max[i];
        if(sstr.fail())
        {
            string_conversion_error_report("mi::math::Bbox<Sint32,3>", value_str);
            return mi::math::Bbox< mi::Sint32, 3 >(0, 0, 0, 0, 0, 0);
        }
    }
    return ret;
}

std::string to_string(mi::math::Bbox< mi::Sint32, 3 > const & value)
{
    std::stringstream sstr;
    sstr << value.min[0] << " " << value.min[1] << " " << value.min[2] << " "
         << value.max[0] << " " << value.max[1] << " " << value.max[2];
    return sstr.str();
}

//----------------------------------------------------------------------
bool load_string_from_file(
    std::string const & infname,
    std::string & data_str,
    std::string & error_mes)
{
    std::string const eh = "Error! _read_file_to_string: ";

    if(infname.empty()){
        error_mes = "load_string_from_file: can not open empty filename file.";
        return false;
    }

    FILE * pfile = fopen(infname.c_str(), "rb");
    if(pfile == 0){
        error_mes = "load_string_from_file: can not open [" + infname + "].";
        return false;
    }

    // get length of file
    fseek(pfile, 0L, SEEK_END);
    mi::Sint64 const file_length = ftell(pfile);
    if(file_length == 0){
        error_mes = "load_string_from_file: empty file.";
        fclose(pfile);
        return false;
    }

    // limit 2G for this function
    if(file_length >= std::numeric_limits< mi::Sint32 >::max()){
        std::stringstream sstr;
        sstr << eh << "too large file size [" << file_length << "] >= "
             << std::numeric_limits< mi::Sint32 >::max();
        error_mes = sstr.str();
        fclose(pfile);
        return false;
    }

    // set buffer size
    data_str.reserve(static_cast< std::string::size_type >(file_length));

    // back to the top of the file
    rewind(pfile);

    // load the file
    while(!feof(pfile)){
        mi::Sint32 const c = fgetc(pfile);
        if(c != EOF){
            data_str += static_cast< char >(c);
        }
    }

    fclose(pfile);

    return true;
}

//----------------------------------------------------------------------
bool get_string_dict_from_inline_parameter_string(
    std::string const & inline_pstring,
    String_dict & out_dict,
    std::string & err_mes
    )
{
    std::string fn = "get_string_dict_from_inline_parameter_string: ";

    std::vector< std::string > key_value_vec;
    std::string const separator = ",";

    nv::index_common::Tokenizer::parse(inline_pstring, separator, key_value_vec);
    String_dict out_dict_tmp;

    // get each key_value
    size_t len = key_value_vec.size();
    for(size_t i = 0; i < len; ++i){
        std::vector< std::string > key_value;
        nv::index_common::Tokenizer::parse(key_value_vec[i], "=", key_value);
        if(key_value.size() != 2){
            // no key=value (or missing one of them)
            std::stringstream sstr;
            sstr << "Error! " << fn << "[" << i << "]-th element ["
                 << key_value_vec[i] << "] is not a 'key=value'.";
            err_mes = sstr.str();
            return false;
        }
        if(out_dict_tmp.is_defined(key_value[0])){
            // duplicated key has no sense for this. Cast an error.
            std::stringstream sstr;
            sstr << "Error! " << fn << "[" << i << "]-th element ["
                 << key_value_vec[i] << "]'s key is duplicated.";
            err_mes = sstr.str();
            return false;
        }
        out_dict_tmp.insert(key_value[0], key_value[1]);
    }

    out_dict = out_dict_tmp;

    return true;
}

//----------------------------------------------------------------------
bool is_valid_key(std::string const & key_str)
{
    const std::string invalid_char = " \t";
    const size_t invalid_char_length = invalid_char.length();
    
    for(size_t i = 0; i < invalid_char_length; ++i){
        mi::Uint8 ch = invalid_char.at(i);
        if(key_str.find(ch) != std::string::npos){
            // found an invalid char
            return false;
        }
    }
    return true;
}

//----------------------------------------------------------------------
nv::index_common::String_dict get_string_dict_from_string(const std::string& str)
{
    nv::index_common::String_dict dict;
    std::stringstream sstr(str);
    dict.read(sstr);

    return dict;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common

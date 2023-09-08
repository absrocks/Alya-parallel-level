/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief string dictionary with IO (map: string -> string)

#ifndef NVIDIA_INDEX_BIN_COMMON_STRING_DICT_H
#define NVIDIA_INDEX_BIN_COMMON_STRING_DICT_H

#include <mi/dice.h>
#include <map>
#include <string>
#include <vector>

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
/// string dictionary (map) with IO.
class String_dict
{
private:
    /// error conversion report switch
    static bool S_is_error_conversion_report_on;

public:
    /// Is error conversion report on?
    /// \return true when error conversion report on
    static bool is_error_conversion_report_on() {
        return S_is_error_conversion_report_on;
    }

    /// set error conversion report on?
    /// \param[in] is_on true when error conversion report on
    static void set_error_conversion_report_on(bool is_on) {
        S_is_error_conversion_report_on = is_on;
    }

private:
    /// Type of the internally used \c hash_map
    typedef std::map< std::string, std::string > String_dict_hash_map;

public:
    /// Value type of the set representation. That is
    /// <tt>std::pair<const key_type, value_type></tt>.
    typedef String_dict_hash_map::value_type     value_type;

    /// Const iterator type for the public API.
    typedef String_dict_hash_map::const_iterator const_iterator;

    /// Forward iterator, partially mutable. The dictionary value can
    /// be changed, but the key cannot be changed.
    typedef String_dict_hash_map::iterator       iterator;

    /// Default constructor initializes the dictionary to the empty set.
    String_dict();

    /// Copy constructor.
    String_dict(String_dict const & dict);

    /// Destructor
    ~String_dict();

    /// Insert an string with the given key and value. If a value of
    /// the same key is already stored in the set, the insertion
    /// replaces the value of that value with the new value.
    ///
    /// \param[in] key    key
    /// \param[in] value value
    /// \return pair of iterator and bool. iterator indicates newly
    /// inserted element, and bool is insertion is succeeded or
    /// not. When the dictionary has already the same key, the bool
    /// returns false. Same as std::hash_map's return value.
    std::pair< String_dict::iterator, bool> insert(std::string const & key,
                                                   std::string const & value);

    /// Insert all entries of \p dict to this.
    /// Each element of \p is inserted.
    /// \return number of newly inserted entries.
    mi::Sint32 insert_all(const String_dict& dict);

    /// Insert the entries of \p dict which are not already defined in
    /// the current dictionary.
    ///
    /// \return number of newly inserted entries.
    mi::Sint32 insert_new(const String_dict& dict);

    /// Erase the dictionary_value with the given \p key. Has no
    /// effect if the \p key does not exist.
    void erase(const std::string & key);

    /// get the value associated with a key. If not found the key,
    /// default_value will be returned.
    ///
    /// \param[in] key key
    /// \param[in] default_value returning default value when key
    /// doesn't exist.
    /// \return value associated with key, when key doen't exist
    /// default_value.
    std::string get(std::string const & key,
                    std::string const & default_value = std::string()) const;

    /// set the value with associated key. The (key, value) exists or
    /// not doesn't matter.
    ///
    /// \param[in] key   key
    /// \param[in] value value
    void set(std::string const & key,
             std::string const & value );

    /// Clear all entries.
    void clear();

    /// Return \c true if the dictionary_value with the given \p key
    /// exists in the set.  Returns always \c false for internal keys.
    bool is_defined(std::string const & key) const;

    /// Return the number of stored elements.
    size_t size() const;

    /// Return a const-iterator pointing to the beginning of the
    /// element set.
    const_iterator begin() const;

    /// Return a const-iterator pointing past-the-end of the
    /// element set.
    const_iterator end() const;

    /// write the entries to \p os.
    /// \param[in] os     output stream
    /// \param[in] prefix prefix of the key
    /// \param[in] is_sort when true, the output is sorted by std::sort.
    void write(std::ostream & os,
               std::string const & prefix = std::string(""),
               bool is_sort = false) const;

    /// read the entries from \p is stream. A line starts with '#'
    /// treated as comment line and ignored.  An empty line (only
    /// contains white spaces) is also ignored.
    ///
    /// \param[in] is input stream. usually stringstream in neuray.
    void read(std::istream & is);

    /// inequality operator
    /// \param[in]  rhs right hand side
    bool operator!=(const String_dict &rhs ) const;

    /// inequality operator
    /// \param[in]  rhs right hand side
    bool operator==(const String_dict &rhs ) const;

private:
    String_dict_hash_map m_map;
};

//----------------------------------------------------------------------
/// check if all keys are defined in the dictionary.
/// string list version.
/// \param[in]  option       checking this option
/// \param[in]  key_list     key list
/// \param[out] p_undef_list (output) when non 0, undefined key list is returned in
/// here.
/// \return true all keys are defined in the option.
bool is_all_keys_defined(
    String_dict const & option,
    std::vector< std::string > const & key_list,
    std::vector< std::string > * p_undef_list = 0
    );

//----------------------------------------------------------------------
/// check if all keys are defined in the dictionary.
/// const char *[] version. The array must 0 terminated.
/// \param[in]  option       checking this option
/// \param[in]  p_key_list   key list. needed 0 terminated.
/// \param[out] p_undef_list (output) when non 0, undefined key list is returned in
/// here.
/// \return true all keys are defined in the option.
bool is_all_keys_defined(
    String_dict const & option,
    char const * const p_key_list[],
    std::vector< std::string > * p_undef_list = 0
    );

//----------------------------------------------------------------------
/// check if all keys are defined in the dictionary and undefined list report by a string.
/// string list version.
/// \param[in]  option       checking this option
/// \param[in]  key_list     key list
/// \param[out] p_undef_list_string (output) when non 0, undefined key
/// list is returned in a string form.
/// \return true all keys are defined in the option.
bool is_all_keys_defined_with_string_report(
    String_dict const & option,
    std::vector< std::string > const & key_list,
    std::string & undef_list_str
    );

//----------------------------------------------------------------------
/// check if all keys are defined in the dictionary and undefined list report by a string.
/// const char *[] version. The array must 0 terminated.
/// \param[in]  option       checking this option
/// \param[in]  p_key_list   key list. needed 0 terminated.
/// \param[out] p_undef_list (output) when non 0, undefined key list is returned in
/// here.
/// \return true all keys are defined in the option.
bool is_all_keys_defined_with_string_report(
    String_dict const & option,
    char const * const p_key_list[],
    std::string & undef_list_str
    );

//----------------------------------------------------------------------
/// argv[] like string array to dictionary converter.
///
/// The string array is the sequence of pair of '-key value' or
/// '-key'.
///
///  - the first string is always has key "command:" .
///  - the pattern is not the -key, its key is "arg_n:".
///  - "command_line:" keeps whole the command line.
///
/// Example:
///  "dict" "-input" "test.mi" "-output" "out.mi" "-verbose" "-opt" "4" "extraarg1"
/// will be:
///
///  - key="command:", value="dict"
///  - key="input",    value="test.mi"
///  - key="output",   value="out.mi"
///  - key="verbose",  value="&lt;no_value&gt;" ... option.is_defined("verbose") is true.
///  - key="opt",      value="4"
///  - key="arg_0:",   value="extraarg1"
///
/// When set_true_if_key_only is true, and only option is defined,
/// then the value becomes true. This is for when preset verbose = 0,
/// then option has -verbose, then the value of "verbose" becomes
/// true.  When set_true_if_key_only is false, if the value is not
/// defined, then "<no_value>" is the value.
///
/// \param[in]  argc argc of main
/// \param[in]  argv argv of main
/// \param[out] option (output) option to be set
/// \param[in]  set_true_if_key_only This is for when the key only defined
/// case. There is no effect if value is found. When true, value is
/// 1. when false, the value is &lt;no_value&gt;.
void string_array_to_string_dict(
    mi::Sint32  const    argc,
    char *        argv[],
    String_dict & option,
    bool set_true_if_key_only = false
    );

//----------------------------------------------------------------------
/// output stream operator for String_dict
std::ostream & operator<<(std::ostream & os,
                          String_dict const & dict);

//----------------------------------------------------------------------
/// magic/header check
///
/// Analyze the first line of the file.
/// Example:
/// "#! String_dict 0"
///
/// \param[in]  fname        filename to check
/// \param[in]  magic_string magic string (Ex. "String_dict")
/// \param[out] file_version (output) version number as a string
/// \param[out] p_error_mes  (output) error message when false
/// returns. If 0, no message stored.
/// \return 0 ... success, -1 ... file not found, -2 ... illegal header
mi::Sint32 get_file_magic(std::string const & fname,
                          std::string const & magic_string,
                          mi::Sint32 &        file_version,
                          std::string *       p_error_mes);

//----------------------------------------------------------------------
/// get string_dict from file with magic check.
///
/// \param[in] dict_opt_fname string dictionary option filename
/// \param[in] magic_str      magic string in the file header
/// \param[out] dict_opt      (output) option
/// \param[out] file_version  (output) file version in the header
/// \param[out] p_error_mes   (output) error message when false
/// returns. If 0, no message stored.
/// \return true when succeeded
bool get_string_dict_with_magic(std::string const & dict_opt_fname,
                                std::string const & magic_str,
                                nv::index_common::String_dict & dict_opt,
                                mi::Sint32        & file_version,
                                std::string       * p_error_mes);

//----------------------------------------------------------------------
/// string_dict key prefix filter.
///
/// This function filters out the entries of the input dictionary
/// according to its keys. This is similer with grep "prefix::" on the
/// key.
///
/// \param[in]  src_dict input source dictionary
/// \param[in]  prefix   prefix string to filter
/// \param[out] filtered filtered result dictionary
/// \param[in]  is_delete_prefix when true, the prefix are removed.
/// \return number of filtered entries. -1 when error.
mi::Sint32 string_dict_key_prefix_filter(
    String_dict const & src_dict,
    std::string const & prefix,
    String_dict       & filtered,
    bool                is_delete_prefix = false);

//----------------------------------------------------------------------
// type conversion helpers

/// get bool from an input string
///
/// \param[in] value_str input value string
/// \return bool value
bool get_bool(std::string const & value_str);

/// to string from bool
///
/// \param[in] value input value
/// \return string representation
std::string to_string(bool value);

//----------------------------------------------------------------------
/// get Sint32 from an input string
/// \param[in] value_str input value string
/// \return Sint32 value
mi::Sint32 get_sint32(std::string const & value_str);

/// to string from mi::Sint32
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Sint32 value);

//----------------------------------------------------------------------
/// get Uint32 from an input string
/// \param[in] value_str input value string
/// \return Uint32 value
mi::Uint32 get_uint32(std::string const & value_str);

/// to string from mi::Uint32
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Uint32 value);

//----------------------------------------------------------------------

/// get Sint64 from an input string
/// \param[in] value_str input value string
/// \return Sint64 value
mi::Sint64 get_sint64(std::string const & value_str);

/// to string from mi::Sint64
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Sint64 value);

//----------------------------------------------------------------------

/// get Uint64 from an input string
/// \param[in] value_str input value string
/// \return Uint64 value
mi::Uint64 get_uint64(std::string const & value_str);

/// to string from mi::Uint64
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Uint64 value);

//----------------------------------------------------------------------

/// get Float32 from an input string
/// \param[in] value_str input value string
/// \return Float32 value
mi::Float32 get_float32(std::string const & value_str);

/// to string from mi::Float32
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Float32 value);

//----------------------------------------------------------------------
/// get Float64 from an input string
/// \param[in] value_str input value string
/// \return Float64 value
mi::Float64 get_float64(std::string const & value_str);

/// to string from mi::Float64
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::Float64 value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Sint32, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Sint32, 2> value
mi::math::Vector_struct< mi::Sint32, 2 > get_vec_sint32_2(std::string const & value_str);

/// to string from mi::math::Vector<mi::Sint32, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Sint32, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Sint64, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Sint64, 2> value
mi::math::Vector_struct< mi::Sint64, 2 > get_vec_sint64_2(std::string const & value_str);

/// to string from mi::math::Vector<mi::Sint64, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Sint64, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Sint32, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Sint32, 3> value
mi::math::Vector_struct< mi::Sint32, 3 > get_vec_sint32_3(std::string const & value_str);

/// to string from mi::math::Vector<mi::Sint32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Sint32, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Uint32, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Uint32, 2> value
mi::math::Vector_struct< mi::Uint32, 2 > get_vec_uint32_2(std::string const & value_str);

/// to string from mi::math::Vector<mi::Uint32, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Uint32, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Uint32, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Uint32, 3> value
mi::math::Vector_struct< mi::Uint32, 3 > get_vec_uint32_3(std::string const & value_str);

/// to string from mi::math::Vector<mi::Uint32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Uint32, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Sint64, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Sint64, 3> value
mi::math::Vector_struct< mi::Sint64, 3 > get_vec_sint64_3(std::string const & value_str);

/// to string from mi::math::Vector<mi::Sint64, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Sint64, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Float32, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Float32, 2> value
mi::math::Vector_struct< mi::Float32, 2 > get_vec_float32_2(std::string const & value_str);

/// to string from mi::math::Vector<mi::Float32, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Float32, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Float32, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Float32, 3> value
mi::math::Vector_struct< mi::Float32, 3 > get_vec_float32_3(std::string const & value_str);

/// to string from mi::math::Vector<mi::Float32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Float32, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Vector_struct<mi::Float32, 4> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Vector_struct<mi::Float32, 4> value
mi::math::Vector_struct< mi::Float32, 4 > get_vec_float32_4(std::string const & value_str);

/// to string from mi::math::Vector<mi::Float32, 4>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Vector< mi::Float32, 4 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Matrix_struct<mi::Float32, 4, 4 > from an input string
/// \param[in] value_str input value string
/// \return mi::math::Matrix_struct<mi::Float32, 4, 4 > value
mi::math::Matrix_struct< mi::Float32, 4, 4 > get_mat_float32_4_4(
    std::string const & value_str);

/// to string from mi::math::Matrix<mi::Float32, 4, 4 >
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Matrix< mi::Float32, 4, 4 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Matrix_struct<mi::Float64, 4, 4 > from an input string
/// \param[in] value_str input value string
/// \return mi::math::Matrix_struct<mi::Float64, 4, 4 > value
mi::math::Matrix_struct< mi::Float64, 4, 4 > get_mat_float64_4_4(
    std::string const & value_str);

/// to string from mi::math::Matrix<mi::Float64, 4, 4 >
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Matrix< mi::Float64, 4, 4 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Color_struct from an input string
/// \param[in] value_str input value string
/// \return mi::math::Color_struct value
mi::math::Color_struct get_color(std::string const & value_str);

/// to string from mi::math::Color
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Color const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Float32, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Float32, 2> value
mi::math::Bbox_struct< mi::Float32, 2 > get_bbox_float32_2(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Float32, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Float32, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Float64, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Float64, 2> value
mi::math::Bbox_struct< mi::Float64, 2 > get_bbox_float64_2(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Float64, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Float64, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Float32, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Float32, 3> value
mi::math::Bbox_struct< mi::Float32, 3 > get_bbox_float32_3(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Float32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Float32, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Uint32, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Uint32, 3> value
mi::math::Bbox_struct< mi::Uint32, 3 > get_bbox_uint32_3(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Uint32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Uint32, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Sint64, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Sint64, 3> value
mi::math::Bbox_struct< mi::Sint64, 3 > get_bbox_sint64_3(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Sint64, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Sint64, 3 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Sint32, 2> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Sint32, 2> value
mi::math::Bbox_struct< mi::Sint32, 2 > get_bbox_sint32_2(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Sint32, 2>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Sint32, 2 > const & value);

//----------------------------------------------------------------------
/// get mi::math::Bbox_struct<mi::Sint64, 3> from an input string
/// \param[in] value_str input value string
/// \return mi::math::Bbox_struct<mi::Sint32, 3> value
mi::math::Bbox_struct< mi::Sint32, 3 > get_bbox_sint32_3(std::string const & value_str);

/// to string from mi::math::Bbox<mi::Sint32, 3>
/// \param[in] value input value
/// \return string representation
std::string to_string(mi::math::Bbox< mi::Sint32, 3 > const & value);

//----------------------------------------------------------------------
/// load string from a file. stdio version (no input file stream version).
///
/// Note: this function keep all the file contents in the memory,
/// therefore it is not good for a large file.
///
/// \param[in]  infname
/// \param[out] data_str
/// \param[out] error_mes error message when return false
/// \return true when file is loaded to data_str.
bool load_string_from_file(
    std::string const & infname,
    std::string & data_str,
    std::string & error_mes);

//----------------------------------------------------------------------
/// get a String dict from a inline parameter string
///
/// Inline parameter string is one line string with comma separated key=value.
/// E.g., "key=value,key=value,..."
///
/// \param[in]  inline_pstring inline parameter string
/// \param[out] out_dict output String_dictionary
/// \param[out] err_mes  error message if returns false
/// \return true when file parsing succeeded
bool get_string_dict_from_inline_parameter_string(
    std::string const & inline_pstring,
    String_dict & out_dict,
    std::string & err_mes
    );

//----------------------------------------------------------------------
/// Test the key validity.
///
/// No white spaces in the key.
///
/// \param[in] key_str key string to be checked.
/// \return true when file parsing succeeded
bool is_valid_key(std::string const & key_str);

//----------------------------------------------------------------------
/// get string_dict from a string
///
/// \param[in] str a string, should be string dict format
/// \return a string dict representation of str.
nv::index_common::String_dict get_string_dict_from_string(const std::string& str);

//----------------------------------------------------------------------
}} // namespace nv::index_common

#endif  // NVIDIA_INDEX_BIN_COMMON_STRING_DICT_H

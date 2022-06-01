/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief string Tokenizer
///
/// A simple string tokenizer. A string is tokenized with separator
/// character set.
///
#ifndef INDEX_BIN_COMMON_TOKENIZER_H
#define INDEX_BIN_COMMON_TOKENIZER_H

#include <string>
#include <vector>

namespace nv {
namespace index_common {

/// Simple string tokenizer class.
/// Examples are in \see test_string_tokenizer.cpp
class Tokenizer {
public:
    /// Constructor.
    /// \param[in] s the string to tokenize
    /// \param[in] separators the separator characters as a pointer to chars
    /// \param[in] return_separators_as_tokens request separators to be returned as tokens
    Tokenizer(
        std::string const& s,
        char const *       separators,
        bool               return_separators_as_tokens = false)
    {
        this->init(s, std::string(separators ? separators : ""), return_separators_as_tokens);
    }

    /// Constructor.
    /// \param[in] s the string to tokenize
    /// \param[in] separators the separator characters as a ref to a string
    /// \param[in] return_separators_as_tokens request separators to be returned as tokens
    Tokenizer(
        std::string const & s,
        std::string const & separators,
        bool                return_separators_as_tokens = false)
    {
        this->init(s, separators, return_separators_as_tokens);
    }

    /// Check if there are more tokens.
    /// \return true if more tokens available
    bool has_tokens() const
    {
        return m_head < m_string.size();
    }

    /// Get the current token.
    /// \return the current token
    std::string token() const
    {
        return m_string.substr(m_head, m_tail - m_head);
    }

    /// Get the next token.
    /// \return next token
    std::string next_token()
    {
        std::string t = token();
        if (m_return_separators) {
            m_head = m_tail;
            m_tail = this->find_bounds();
            if (m_head == m_tail)
                m_tail++;
        } else {
            m_head = m_tail + 1;
            m_tail = this->find_bounds();
        }
        return t;
    }

    /// Parse and get a token list.
    /// A convenient function for getting a token list.
    ///
    /// Example:
    /// \code
    ///   std::vector< std::string > token_list;
    ///   Tokenizer::parse("A Non-Manifold,Mesh", " ,", token_list);
    /// \endcode
    /// gives
    ///   token_list[0] == "A";
    ///   token_list[1] == "Non-Manifold";
    ///   token_list[2] == "Mesh";
    ///
    /// Notice: If the last entry is empty, it does not count.
    ///
    /// For example, tokenize "hello," with "," gives only one token {
    /// "hello" } instead of two tokens {"hello", ""}. But, tokenize
    /// "hello,,world" with "," gives three tokens { "hello", "",
    /// "world" }.
    ///
    /// \param[in] source_str the source string
    /// \param[in] separators the separators
    /// \param[in,out] token_list the tokens will end up here
    static void parse(
        std::string const & source_str,
        std::string const & separators,
        std::vector< std::string >& token_list)
    {
        Tokenizer st(source_str, separators);
        while(st.has_tokens()){
            token_list.push_back(st.token());
            st.next_token();
        }
    }

    /// find_first_of, but with checking 'what' is really in there.
    /// \param[in] source source string to be searched
    /// \param[in] what   string we are looking for
    /// \param[in] offset current looking at position
    static std::string::size_type tokenizer_real_find_first_of(
        std::string const & source,
        std::string const & what,
        std::string::size_type offset)
    {
        std::string::size_type pos = source.find_first_of(what, offset);
        if (pos != std::string::npos) {
            // check that the complete string what fits into
            if (pos + what.size() >= source.size()) {
                pos = std::string::npos;
            }
            else {
                // check that we have found a real occurrence of what
                for (size_t i=0; i<what.size(); ++i) {
                    if (what[i] == source[pos+i]){
                        continue;
                    }
                    pos = tokenizer_real_find_first_of(source, what, offset+1);
                }
            }
        }
        return pos;
    }

    /// Retrieve the tokenized string when dividing the \p source at the \p separator. Note that
    /// this version matches the \b complete separator, not single characters of it. I.e.
    /// \code
    ///  string source("A::B:::C");
    ///  vector<string> res = Tokenizer::tokenize_string(source, "::");
    ///  assert(res[0] == "A" && res[1] == "B" && res[2] == ":C");
    /// \endcode
    /// \param[in] source     the input string
    /// \param[in] separator  separator string
    /// \return list of divided input
    static std::vector<std::string> tokenize_string(
        std::string const & source,
        std::string const & separator)
    {
        std::vector<std::string> result;

        std::string::size_type pos = 0;
        std::string::size_type last = pos;
        while (pos < source.size()
               && (pos = tokenizer_real_find_first_of(source, separator, last)) != std::string::npos)
        {
            result.push_back(source.substr(last, pos-last));
            pos += separator.size();
            last = pos;
        }

        if (last < source.size()){
            result.push_back(source.substr(last));
        }
        return result;
    }

private:
    /// Common init for constructors
    /// \param[in] s the string to tokenize
    /// \param[in] separators the separator characters
    /// \param[in] return_separators_as_tokens request separators to be returned as tokens
    void init(
        std::string const & s,
        std::string const & separators,
        bool                return_separators_as_tokens)
     {
         m_string = s;
         m_separators = separators;
         m_head = 0;
         m_tail = this->find_bounds();         
         m_return_separators = return_separators_as_tokens;
     }

    /// Get the index of the next token.
    /// \return index of next token
    size_t find_bounds() const
    {
        size_t tail = m_string.size();
        size_t num_separators = m_separators.size();
        for (size_t i = m_head; i < tail; i++){
            for (size_t j = 0; j < num_separators; j++){
                if (m_string[i] == m_separators[j]){
                    return i;
                }
            }
        }
        return tail;
    }

private:
    std::string     m_string;               ///< the string to tokenize
    std::string     m_separators;           ///< the string of separators
    size_t          m_head;                 ///< start index of current token
    size_t          m_tail;                 ///< end index of current token
    bool            m_return_separators;    ///< return separators as tokens
};


}} // namespace nv::index_common

#endif // INDEX_BIN_COMMON_TOKENIZER_H

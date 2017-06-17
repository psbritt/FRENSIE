//---------------------------------------------------------------------------//
//!
//! \file   Utility_FromStringTraits.hpp
//! \author Alex Robinson
//! \brief  FromStringTraits class specializations
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_FROM_STRING_TRAITS_HPP
#define UTILITY_FROM_STRING_TRAITS_HPP

// Std Lib Includes
#include <sstream>
#include <string>
#include <type_traits>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>
#include <list>
#include <forward_list>
#include <deque>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

// Boost Includes
#include <boost/algorithm/string/trim.hpp>

// FRENSIE Includes
#include "Utility_ExceptionTestMacros.hpp"

namespace Utility{
  
// Initialize the input stream that will be used to deserialize an object
inline void initializeInputStream( std::istream& is, const char start_delim );
{
  // Read the start deliminator
  std::string start_bracket;
  std::getline( is, start_bracket, start_delim );
  start_bracket = Teuchos::Utils::trimWhiteSpace( start_bracket );

  TEST_FOR_EXCEPTION( start_bracket.size() != 0,
                      std::runtime_error,
                      "Unable to initialize the input stream because there "
                      "are characters preceding the start deliminator ("
                      << start_delim << ")! Here are the extra characters: "
                      << start_bracket );
}

// Move the input stream to the start of the next element
inline bool moveInputStreamToNextElement( std::istream& is,
                                          const char elem_delim,
                                          const char end_delim )
{
  // Search for the specified element deliminator
  char delim;
  
  while( true )
  {
    is.get( delim );

    TEST_FOR_EXCEPTION( !is,
                        std::runtime_error,
                        "Unable to move the input stream to the next element "
                        "(EOF reached unexpectedly)!" );

    // Another element must be deserialized
    if( delim == elem_delim )
      return false;
    
    // All elements have been deserialized
    else if( delim == end_delim )
      return true;

    // An invalid deliminator has been encountered
    else if( delim != ' ' )
    {
      THROW_EXCEPTION( std::runtime_error,
                       "Unable to move the input stream to the next element "
                       "(bad deliminator encountered: " << delim << ")!" );
      return true;
    }
  }
}

/*! Specialization of FromStringTraits for std::string
 * \ingroup from_string_traits
 */
template<>
struct FromStringTraits<std::string>
{
  //! Return the string
  static inline std::string fromString( const std::sting& obj_rep )
  { return obj_rep; }

  //! Extract the string from the stream
  static inline void fromStream( std::istream& is,
                                 std::string& obj,
                                 const std::string& delims = std::string() )
  {
    std::string obj;

    if( delims.size() > 0 )
    {
      bool done = false;
      
      while( !done )
      {
        char string_element;
        is.get( string_element );

        TEST_FOR_EXCEPTION( !is,
                            std::runtime_error,
                            "Unable to get the string element (EOF reached "
                            "unexpectedly)!" );

        for( size_t i = 0; i < delims.size(); ++i )
        {
          // A deliminator has been reached - stop reading from stream
          if( string_element == delims[i] )
          {
            done = true;

            // Put the deliminator back in the stream so that it can be
            // parsed correctly later
            is.putback( string_element );
            
            break;
          }
        }

        if( !done )
          obj.push_back( string_element );
      }
    }
    // No deliminators have been specified - use a space
    else
      std::getline( is, obj, ' ' );

    // Trim the extracted string
    boost::algorithm::trim( obj );
  }
};

/*! Specialization of FromStringTraits for bool
 * \ingroup from_string_traits
 */
template<>
struct FromStringTraits<bool>
{
  //! Convert the string to a boolean
  static inline bool fromString( const std::string& obj_rep )
  {
    std::string lower_case_obj_rep =
      boost::algorithm::to_lower_copy( obj_rep );

    boost::algorithm::trim( lower_case_obj_rep );

    if( lower_case_obj_rep == "true" || lower_case_obj_rep == "1" )
      return true;
    else if( lower_case_obj_rep == "false" || lower_case_obj_rep == "0" )
      return false;
    else
    {
      THROW_EXCEPTION( std::runtime_error,
                       "The string (" << obj_rep << ") does not correspond to "
                       "a boolean!" );
      return false;
    }
  }

  //! Extract a boolean from the stream
  static inline void fromStream( std::istream& is,
                                 bool& obj,
                                 const std::string& delim = std::string() )
  {
    std::string bool_rep;

    Utility::fromStream( is, bool_rep, delim );

    return FromStringTraits<bool>::fromString( bool_rep );
  }
};

/*! Partial specialization of FromStringTraits for arithmetic types
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<T,typename std::enable_if<std::is_arithmetic<T>::value>::type>
{
  //! Convert the string to an object of type T
  static inline T fromString( const std::string& obj_rep )
  {
    std::istringstream iss( obj_rep );

    T obj;

    FromStringTraits<T>::fromStream( is, obj );

    return obj;
  }

  //! Extract the object from a stream
  static inline void fromStream( std::istream& is,
                                 T& obj,
                                 const std::string& = std::string() )
  { is >> obj; }
};

namespace Details{

// The helper class for stream i/o ops with tuples
template<size_t I, typename TupleType, typename Enable = void>
struct TupleFromStreamHelper
{
  static inline void fromStream( std::istream& is, TupleType& tuple )
  {
    // Extract the tuple element
    try{
      Utility::fromStream( is, std::get<I>( tuple ), ",}" );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Tuple element " << I << " was not successfully "
                             "extracted from the input stream!" );

    // Position the stream at the start of the next element (or end)
    try{
      Utility::moveInputStreamToNextElement( is, ',', '}' );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Could not move the input stream to the next "
                             "element (last tuple element successfully "
                             "extracted = " << I << ")!" );

    // Extract the remaining tuple elements
    TupleFromStreamHelper<I+1,TupleType>::fromStream( is, tuple );
  }
};

/*! \brief The helper class for stream i/o ops with tuples
 * (specialization for I == std::tuple_size-1, which is for the last element)
 */
template<size_t I, typename TupleType>
struct TupleFromStreamHelper<I, TupleType, typename std::enable_if<I==std::tuple_size<TupleType>::value-1>::type>
{
  static inline void fromStream( std::istream& is, TupleType& tuple )
  {
    // Extract the tuple element
    try{
      Utility::fromStream( is, std::get<I>( tuple ), ",}" );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Tuple element " << I << " was not "
                             "successfully extracted from the input stream!" );

    // Position the stream at the start of the next element (or end)
    bool at_end;
    try{
      at_end = Utility::moveInputStreamToNextElement( is, ',', '}' );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Could not move the input stream to the next "
                             "element (last tuple element successfully "
                             "extracted = " << I << ")!" );

    TEST_FOR_EXCEPTION( !at_end,
                        std::runtime_error,
                        "Finished extracting tuple before the end of the "
                        "stream was reached (this likely means that the "
                        "string is not compatible with the tuple type)!" );
  }
};

/*! \brief The helper class for stream i/o ops with tuples
 * (specialization for I == std::tuple_size, which is past the last element)
 */
template<size_t I, typename TupleType>
struct TupleFromStreamHelper<I, TupleType, typename std::enable_if<I==std::tuple_size<TupleType>::value>::type>
{
  static inline void fromStream( std::istream& is, TupleType& tuple )
  { /* ... */ }
};

} // end Details namespace

/*! Partial specialization of FromStringTraits for std::tuple
 * \ingroup from_string_traits
 */
template<typename... Types>
struct FromStringTraits<std::tuple<Types...> >
{
  //! Convert the string to an object of type T
  static inline std::tuple<Types...> fromString( const std::string& obj_rep )
  {
    std::istringstream iss( obj_rep );

    std::tuple<Types...> obj;

    FromStringTraits<std::tuple<Types...> >::fromStream( is, obj );

    return obj;
  }

  //! Extract the object from a stream
  static inline void fromStream( std::istream& is,
                                 std::tuple<Types...>& obj,
                                 const std::string& = std::string() )
  { 
    try{
      // Initialize the input stream
      Utility::initializeInputStream( is, '{' );

      Details::TupleFromStreamHelper<0,std::tuple<Types...> >::fromStream( is, obj );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Could not extract a tuple from the stream!" );
  }
};

/*! Partial specialization of FromStringTraits for std::pair
 * \ingroup from_string_traits
 */
template<typename T1, typename T2>
struct FromStringTraits<std::pair<T1,T2> >
{
  //! Convert the string to an object of type T
  static inline std::pair<T1,T2> fromString( const std::string& obj_rep )
  {
    std::pair<T1,T2> pair;

    std::istringstream iss( obj_rep );

    FromStringTraits<std::pair<T1,T2> >::fromStream( is, pair );

    return pair;
  }

  //! Extract the object from a stream
  static inline void fromStream( std::istream& is,
                                 std::pair<T1,T2>& obj,
                                 const std::string& = std::string() )
  { 
    auto pair_reference = std::tie( obj.first, obj.second );

    Utility::fromStream( is, pair_reference );
  }
};

namespace Details{

//! FromStringTraits base class for stl compliant containers
template<typename FromStringTraitsChild>
struct FromStringTraitsSTLCompliantContainerBase
{
  //! Convert the string to the required contianer type
  static inline typename FromStringTraitsChild::ContainerType fromString(
                                                   const std::string& obj_rep )
  {
    typename FromStringTraitsChild::ContainerType container;

    std::istringstream iss( obj_rep );

    FromStringTraitsChild::fromStream( iss, container );

    return container;
  }

protected:

  //! Extract the object from a stream
  template<ElementInsertionMemberFunction>
  static inline void fromStreamImpl(
                            std::istream& is,
                            typename FromStringTraitsChild::ContainerType& obj,
                            ElementInsertionMemberFunction& insert_element )
  {
    try{
      // Initialize the input stream
      Utility::initializeInputStream( is, '{' );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Could not extract a tuple from the stream!" );

    // Extract each element of the array
    bool done = false;
    size_t element_index = 0;

    while( !done )
    {
      typename FromStringTraitsChild::ContainerType::value_type element;

      try{
        element = Utility::fromStream( is, element, ",}" );
      }
      EXCEPTION_CATCH_RETHROW( std::runtime_error,
                               "Element " << element_index << " was not "
                               "successfully extracted from the stream!" );

      // Check if the stream is still valid
      TEST_FOR_EXCEPTION( !is,
                          std::runtime_error,
                          "An error has occurred in the stream while "
                          "attempting to extract element "
                          << element_index << "!" );

      // The element has been successfully extracted
      insert_element( obj, element );
      ++element_index;

      // Position the stream at the start of the next element (or end)
      try{
        done = Details::moveInputStreamToNextElement( is, ',', '}' );
      }
      EXCEPTION_CATCH_RETHROW( std::runtime_error,
                               "Could not move the input stream to the next "
                               "element (" << element_index << ")!" );
    }
  }
};

} // end Details namespace

/*! Partial specialization of FromStringTraits for std::vector
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::vector<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::vector<T> >
{
  //! The container type (used by the base class)
  typedef std::vector<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::push_back );
  }
};

/*! Partial specialization of FromStringTraits for std::list
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::list<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::list<T> >
{ 
  //! The container type (used by the base class)
  typedef std::list<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::push_back );
  }
};

/*! Partial specialization of FromStringTraits for std::forward_list
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::forward_list<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::forward_list<T> >
{ 
  //! The container type (used by the base class)
  typedef std::forward_list<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    // Elements must be appended to the front of a forward list. To preserve
    // the element order in the stream we will have to call reverse once
    // all elements have been extracted
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::push_front );

    obj.reverse();
  }
};

/*! Partial specialization of FromStringTraits for std::deque
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::deque<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::deque<T> >
{ 
  //! The container type (used by the base class)
  typedef std::deque<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::push_back );
  }
};

/*! Partial specialization of FromStringTraits for std::set
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::set<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::set<T> >
{ 
  //! The container type (used by the base class)
  typedef std::set<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::insert );
  }
};

/*! Partial specialization of FromStringTraits for std::unordered_set
 * \ingroup from_string_traits
 */
template<typename T>
struct FromStringTraits<std::unordered_set<T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::unordered_set<T> >
{ 
  //! The container type (used by the base class)
  typedef std::unordered_set<T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::insert );
  }
};

/*! Partial specialization of FromStringTraits for std::map
 * \ingroup from_string_traits
 */
template<typename Key, typename T>
struct FromStringTraits<std::map<Key,T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::map<Key,T> >
{ 
  //! The container type (used by the base class)
  typedef std::map<Key,T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::insert );
  }
};

/*! Partial specialization of FromStringTraits for std::unordered_map
 * \ingroup from_string_traits
 */
template<typename Key, typename T>
struct FromStringTraits<std::unordered_map<Key,T> > : public Details::FromStringTraitsSTLCompliantContainerBase<FromStringTraits<std::unordered_map<Key,T> >
{ 
  //! The container type (used by the base class)
  typedef std::unordered_map<Key,T> ContainerType;

  //! Extract the container from a stream
  static inline void fromStream( std::istream& is,
                                 ContainerType& obj,
                                 const std::string& = std::string() )
  {
    FromStringTraits<ContainerType>::fromStreamImpl( is, obj, &ContainerType::insert );
  }
};

// Convert the string to an object of type T
template<typename T>
inline void fromString( const std::string& obj_rep, T& obj )
{
  obj = Utility::FromStringTraits<typename std::remove_reference<T>::type>::fromString( obj_rep );
}

// Extract an object of type T from the stream
template<typename T>
inline void fromStream( std::istream& is, T& obj, const std::string& delims )
{
  Utility::FromStringTraits<typename std::remove_reference<T>::type>::fromStream( is, obj, delims );
}
  
} // end Utility namespace

#endif // end UTILITY_FROM_STRING_TRAITS_HPP

//---------------------------------------------------------------------------//
// end Utility_FromStringTraits.hpp
//---------------------------------------------------------------------------//

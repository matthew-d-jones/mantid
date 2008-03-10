//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidKernel/PropertyManager.h"
#include "MantidKernel/Exception.h"
#include <algorithm>

namespace Mantid
{
namespace Kernel
{

// Get a reference to the logger
Logger& PropertyManager::g_log = Logger::get("PropertyManager");

/// Default constructor
PropertyManager::PropertyManager() :
  m_properties(),
  m_orderedProperties()
{
}

/// Virtual destructor
PropertyManager::~PropertyManager()
{
  for ( PropertyMap::iterator it = m_properties.begin(); it != m_properties.end(); ++it )
  {
    delete it->second;
  }
}

/** Add a property to the list of managed properties
 *  @param p The property object to add
 *  @throw Exception::ExistsError if a property with the given name already exists
 *  @throw std::invalid_argument  if the property declared has an empty name.
 */
void PropertyManager::declareProperty( Property *p )
{
  // Get the name of the property and don't permit empty names
  std::string key = p->name();
  if (key.empty()) 
  {
//    delete p;
    throw std::invalid_argument("An empty property name is not permitted");
  }
  
  std::transform(key.begin(), key.end(), key.begin(), toupper);
  if ( m_properties.insert(PropertyMap::value_type(key, p)).second)
  {
    m_orderedProperties.push_back(p);
  }
  else
  {
//    delete p;
    throw Exception::ExistsError("Property with given name already exists", p->name() );
  }
}

/** Specialised version of declareProperty template method to prevent the creation of a
 *  PropertyWithValue of type const char* if an argument in quotes is passed (it will be
 *  converted to a string). The validator, if provided, needs to be a string validator.
   *  @param name The name to assign to the property
   *  @param value The initial value to assign to the property
   *  @param validator Pointer to the (optional) validator. Ownership will be taken over.
   *  @param doc The (optional) documentation string
   *  @throw Exception::ExistsError if a property with the given name already exists
   *  @throw std::invalid_argument  if the name argument is empty
 */
void PropertyManager::declareProperty( const std::string &name, const char* value,
                                       IValidator<std::string> *validator, const std::string &doc )
{
  // Simply call templated method, converting character array to a string
  declareProperty(name, std::string(value), validator, doc);
}

/** Set the ordered list of properties by one string of values.
 *  @param values The list of property values
 *  @throws Exception::NotImplementedError because it isn't, yet
 */
// Care will certainly be required in the calling of this function or it could all go horribly wrong!
void PropertyManager::setProperties( const std::string &values )
{
  throw Exception::NotImplementedError("Coming to an iteration near you soon...");
}

/** Set the value of a property by string
 *  @param name The name of the property (case insensitive)
 *  @param value The value to assign to the property
 *  @throw Exception::NotFoundError if the named property is unknown
 */
void PropertyManager::setPropertyValue( const std::string &name, const std::string &value )
{
  Property *p = getPointerToProperty(name);   // throws NotFoundError if property not in vector
  bool success = p->setValue(value);
  if ( !success ) throw std::invalid_argument("Invalid value for this property");
}

/** Checks whether the named property is already in the list of managed property.
 *  @param name The name of the property (case insensitive)
 *  @return True if the property is already stored
 */
bool PropertyManager::existsProperty( const std::string& name ) const
{
  try 
  {
    getPointerToProperty(name);
    return true;
  }
  catch (Exception::NotFoundError e)
  {
    return false;
  }
}

/** Validates all the properties in the collection
 *  @return True if all properties have a valid value
 */
bool PropertyManager::validateProperties() const
{
  bool allValid = true;
  for ( PropertyMap::const_iterator it = m_properties.begin(); it != m_properties.end(); ++it )
  {
    if ( ! (it->second->isValid()) )
    {
      g_log.error() << "Property \"" << it->first << "\" is not set to a valid value." << std::endl;
      allValid=false;
    }
  }
  return allValid; 
}

/** Get the value of a property as a string
 *  @param name The name of the property (case insensitive)
 *  @return The value of the named property
 *  @throw Exception::NotFoundError if the named property is unknown
 */
std::string PropertyManager::getPropertyValue( const std::string &name ) const
{
  Property *p = getPointerToProperty(name);   // throws NotFoundError if property not in vector
  return p->value();
}

/** Get a property by name
 *  @param name The name of the property (case insensitive)
 *  @return A pointer to the named property
 *  @throw Exception::NotFoundError if the named property is unknown
 */
Property* PropertyManager::getPointerToProperty( const std::string &name ) const
{
  std::string ucName = name;
  std::transform(ucName.begin(), ucName.end(), ucName.begin(), toupper);
  PropertyMap::const_iterator it = m_properties.find(ucName);
  if (it != m_properties.end())
  {
    return it->second;
  }
  throw Exception::NotFoundError("Unknown property", name);
}

/** Get the list of managed properties.
 *  The properties will be stored in the order that they were declared.
 *  @return A vector holding pointers to the list of properties
 */
const std::vector< Property* >& PropertyManager::getProperties() const
{
  return m_orderedProperties;
}

/// @cond

// getValue template specialisations (there is no generic implementation)
// Note that other implementations can be found in Workspace.cpp & Workspace1D/2D.cpp (to satisfy
// package dependency rules).

template<> DLLExport
int PropertyManager::getValue<int>(const std::string &name) const
{
  PropertyWithValue<int> *prop = dynamic_cast<PropertyWithValue<int>*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template<> DLLExport
bool PropertyManager::getValue<bool>(const std::string &name) const
{
  PropertyWithValue<bool> *prop = dynamic_cast<PropertyWithValue<bool>*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template<> DLLExport
double PropertyManager::getValue<double>(const std::string &name) const
{
  PropertyWithValue<double> *prop = dynamic_cast<PropertyWithValue<double>*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template<> DLLExport
std::vector<int> PropertyManager::getValue<std::vector<int> >(const std::string &name) const
{
  PropertyWithValue<std::vector<int> > *prop = dynamic_cast<PropertyWithValue<std::vector<int> >*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template<> DLLExport
std::vector<double> PropertyManager::getValue<std::vector<double> >(const std::string &name) const
{
  PropertyWithValue<std::vector<double> > *prop = dynamic_cast<PropertyWithValue<std::vector<double> >*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template<> DLLExport
std::vector<std::string> PropertyManager::getValue<std::vector<std::string> >(const std::string &name) const
{
  PropertyWithValue<std::vector<std::string> > *prop = dynamic_cast<PropertyWithValue<std::vector<std::string> >*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    throw std::runtime_error("Attempt to assign property of incorrect type");
  }
}

template <> DLLExport
const char* PropertyManager::getValue<const char*>(const std::string &name) const
{
  return getPropertyValue(name).c_str();
}

// This template implementation has been left in because although you can't assign to an existing string
// via the getProperty() method, you can construct a local variable by saying, 
// e.g.: std::string s = getProperty("myProperty")
template <> DLLExport
std::string PropertyManager::getValue<std::string>(const std::string &name) const
{
  return getPropertyValue(name);
}

template <> DLLExport
Property* PropertyManager::getValue<Property*>(const std::string &name) const
{
  return getPointerToProperty(name);
}

// If a string is given in the argument, we can be more flexible
template <>
void PropertyManager::setProperty<std::string>(const std::string &name, const std::string value)
{
  this->setPropertyValue(name, value);
}
/// @endcond

/** Get the value of a property. Allows you to assign directly to a variable of the property's type
 *  (if a supported type).
 *  
 *  *** This method does NOT work for assigning to an existing std::string.
 *      In this case you have to use getPropertyValue() instead.
 *      Note that you can, though, construct a local string variable by writing,
 *      e.g. std::string s = getProperty("myProperty"). ***
 *  
 *  @param name The name of the property
 *  @return The value of the property. Will be cast to the desired type (if a supported type).
 *  @throw std::runtime_error If an attempt is made to assign a property to a different type
 *  @throw Exception::NotFoundError If the property requested does not exist
 */
PropertyManager::TypedValue PropertyManager::getProperty( const std::string &name ) const
{
  return TypedValue(*this, name);
}

} // namespace Kernel
} // namespace Mantid

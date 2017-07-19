//---------------------------------------------------------------------------//
//!
//! \file   Utility_AnalogElasticDistribution.hpp
//! \author Luke Kersting
//! \brief  Analog elastic distribution class declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_ANALOG_ELASTIC_DISTRIBUTION_HPP
#define UTILITY_ANALOG_ELASTIC_DISTRIBUTION_HPP

// Trilinos Includes
#include <Teuchos_Array.hpp>

// FRENSIE Includes
#include "Utility_TabularOneDDistribution.hpp"
#include "Utility_ParameterListCompatibleObject.hpp"
#include "Utility_InterpolationPolicy.hpp"
#include "Utility_Tuple.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_UnitTraits.hpp"

namespace Utility{

/*! The interpolated distribution class declaration
 * \ingroup one_d_distributions
 */
template<typename InterpolationPolicy,
         typename IndependentUnit,
         typename DependentUnit>
class UnitAwareAnalogElasticDistribution : public UnitAwareTabularOneDDistribution<IndependentUnit,DependentUnit>,
                                     public ParameterListCompatibleObject<UnitAwareAnalogElasticDistribution<InterpolationPolicy,IndependentUnit,DependentUnit> >
{

private:

  // The unnormalized cdf quantity
  typedef typename QuantityTraits<typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::DistNormQuantity>::template GetQuantityToPowerType<-1>::type UnnormCDFQuantity;

  // The slope unit traits
  typedef UnitTraits<typename UnitTraits<DependentUnit>::template GetMultipliedUnitType<typename UnitTraits<IndependentUnit>::InverseUnit>::type> SlopeUnitTraits;

  // The slope quantity
  typedef typename SlopeUnitTraits::template GetQuantityType<double>::type SlopeQuantity;

  // The distribution normalization quantity type
  typedef typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::DistNormQuantity DistNormQuantity;

  // Typedef for QuantityTraits<double>
  typedef QuantityTraits<double> QT;

  // Typedef for QuantityTraits<IndepQuantity>
  typedef QuantityTraits<typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::IndepQuantity> IQT;

  // Typedef for QuantityTraits<InverseIndepQuantity>
  typedef QuantityTraits<typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::InverseIndepQuantity> IIQT;

  // Typedef for QuantityTraits<DepQuantity>
  typedef QuantityTraits<typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::DepQuantity> DQT;

  // Typedef for QuantityTraits<DistNormQuantity>
  typedef QuantityTraits<DistNormQuantity> DNQT;

  // Typedef for QuantityTraits<UnnormCDFQuantity>
  typedef QuantityTraits<UnnormCDFQuantity> UCQT;

public:

  //! This distribution type
  typedef UnitAwareAnalogElasticDistribution<InterpolationPolicy,IndependentUnit,DependentUnit> ThisType;

  //! The independent quantity type
  typedef typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::IndepQuantity IndepQuantity;

  //! The inverse independent quantity type
  typedef typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::InverseIndepQuantity InverseIndepQuantity;

  //! The dependent quantity type
  typedef typename UnitAwareOneDDistribution<IndependentUnit,DependentUnit>::DepQuantity DepQuantity;

  //! Default constructor
  UnitAwareAnalogElasticDistribution();

  //! Basic constructor (potentially dangerous)
  UnitAwareAnalogElasticDistribution(
                        const Teuchos::Array<double>& independent_values,
                        const Teuchos::Array<double>& dependent_values,
                        const unsigned& atomic_number );

  //! Constructor
  template<typename InputIndepQuantity, typename InputDepQuantity>
  UnitAwareAnalogElasticDistribution(
                  const Teuchos::Array<InputIndepQuantity>& independent_values,
                  const Teuchos::Array<InputDepQuantity>& dependent_values,
                  const unsigned& atomic_number );

  //! Copy constructor
  template<typename InputIndepUnit, typename InputDepUnit>
  UnitAwareAnalogElasticDistribution( const UnitAwareAnalogElasticDistribution<InterpolationPolicy,InputIndepUnit,InputDepUnit>& dist_instance );

  //! Construct distribution from a unitless dist. (potentially dangerous)
  static UnitAwareAnalogElasticDistribution fromUnitlessDistribution( const UnitAwareAnalogElasticDistribution<InterpolationPolicy,void,void>& unitless_distribution );

  //! Assignment operator
  UnitAwareAnalogElasticDistribution& operator=(
                           const UnitAwareAnalogElasticDistribution& dist_instance );

  //! Destructor
  ~UnitAwareAnalogElasticDistribution()
  { /* ... */ }

  //! Evaluate the distribution
  DepQuantity evaluate( const IndepQuantity indep_var_value ) const;

  //! Evaluate the PDF
  InverseIndepQuantity evaluatePDF( const IndepQuantity indep_var_value ) const;

  //! Evaluate the CDF
  double evaluateCDF( const IndepQuantity indep_var_value ) const;

  //! Return a random sample from the distribution
  IndepQuantity sample() const;

  //! Return a random sample and record the number of trials
  IndepQuantity sampleAndRecordTrials( unsigned& trials ) const;

  //! Return a random sample and bin index from the distribution
  IndepQuantity sampleAndRecordBinIndex( unsigned& sampled_bin_index ) const;

  //! Return a random sample from the distribution at the given CDF value
  IndepQuantity sampleWithRandomNumber( const double random_number ) const;

  //! Return a random sample from the distribution in a subrange
  IndepQuantity sampleInSubrange( const IndepQuantity max_indep_var ) const;

  //! Return a random sample from the distribution at the given CDF value in a subrange
  IndepQuantity sampleWithRandomNumberInSubrange(
                                     const double random_number,
                                     const IndepQuantity max_indep_var ) const;

  //! Return the upper bound of the distribution independent variable
  IndepQuantity getUpperBoundOfIndepVar() const;

  //! Return the lower bound of the distribution independent variable
  IndepQuantity getLowerBoundOfIndepVar() const;

  //! Return the distribution type
  OneDDistributionType getDistributionType() const;

  //! Test if the distribution is continuous
  bool isContinuous() const;

  //! Method for placing the object in an output stream
  void toStream( std::ostream& os ) const;

  //! Method for initializing the object from an input stream
  void fromStream( std::istream& is );

  //! Method for testing if two objects are equivalent
  bool isEqual( const UnitAwareAnalogElasticDistribution& other ) const;

protected:

  //! Copy constructor (copying from unitless distribution only)
  UnitAwareAnalogElasticDistribution( const UnitAwareAnalogElasticDistribution<InterpolationPolicy,void,void>& unitless_dist_instance, int );

  //! Test if the dependent variable can be zero within the indep bounds
  bool canDepVarBeZeroInIndepBounds() const;

  //! Test if the independent variable is compatible with Lin processing
  bool isIndepVarCompatibleWithProcessingType(
                                        const LinIndepVarProcessingTag ) const;
  
  //! Test if the independent variable is compatible with Log processing
  bool isIndepVarCompatibleWithProcessingType(
                                        const LogIndepVarProcessingTag ) const;

  //! Test if the dependent variable is compatible with Lin processing
  bool isDepVarCompatibleWithProcessingType(
                                          const LinDepVarProcessingTag ) const;

  //! Test if the dependent variable is compatible with Log processing
  bool isDepVarCompatibleWithProcessingType(
                                          const LogDepVarProcessingTag ) const;

  //! Set the norm constant
  void setNormConstant( const DistNormQuantity norm_constant );

  //! Set the max CDF value
  void setMaxCDF( const UnnormCDFQuantity max_cdf );

private:

  // Initialize the distribution
  void initializeDistributionFromRawData(
                              const Teuchos::Array<double>& independent_values,
                              const Teuchos::Array<double>& dependent_values );

  // Initialize the distribution
  template<typename InputIndepQuantity, typename InputDepQuantity>
  void initializeDistribution(
                  const Teuchos::Array<InputIndepQuantity>& independent_values,
                  const Teuchos::Array<InputDepQuantity>& dependent_values );

  // Reconstruct original distribution
  void reconstructOriginalDistribution(
                         Teuchos::Array<IndepQuantity>& independent_values,
                         Teuchos::Array<DepQuantity>& dependent_values ) const;

  // Reconstruct original distribution w/o units
  void reconstructOriginalUnitlessDistribution(
                              Teuchos::Array<double>& independent_values,
                              Teuchos::Array<double>& dependent_values ) const;

  // Convert the unitless values to the correct units
  template<typename Quantity>
  static void convertUnitlessValues(
                                 const Teuchos::Array<double>& unitless_values,
                                 Teuchos::Array<Quantity>& quantitites );

  // Return a random sample using the random number and record the bin index
  IndepQuantity sampleImplementation( double random_number,
                                      unsigned& sampled_bin_index ) const;

  // All possible instantiations are friends
  template<typename FriendInterpolationPolicy,
           typename FriendIndepUnit,
           typename FriendDepUnit>
  friend class UnitAwareAnalogElasticDistribution;

  // The distribution type
  static const OneDDistributionType distribution_type = ANALOG_ELASTIC_DISTRIBUTION;

  // The distribution (first = indep_var, second = cdf, third = pdf,
  // fourth = pdf slope): both the pdf and cdf are left unnormalized to
  // prevent altering the grid with log interpolation
  typedef Teuchos::Array<Quad<IndepQuantity,UnnormCDFQuantity,DepQuantity,SlopeQuantity> > DistributionArray;
  DistributionArray d_distribution;

  // The normalization constant
  DistNormQuantity d_norm_constant;

  // The max CDF value
  UnnormCDFQuantity d_max_cdf;
};

/*! The analog elastic distribution (unit-agnostic)
 * \ingroup one_d_distributions
 */
template<typename InterpolationPolicy> using AnalogElasticDistribution =
  UnitAwareAnalogElasticDistribution<InterpolationPolicy,void,void>;

} // end Utility namespace

namespace Teuchos{

/*! Type name traits specialization for the Utility::AnalogElasticDistribution
 *
 * \details The name function will set the type name that must be used in
 * xml files.
 */
template<typename InterpolationPolicy>
class TypeNameTraits<Utility::AnalogElasticDistribution<InterpolationPolicy> >
{
public:
  static std::string name()
  {
    std::ostringstream iss;
    iss << "Analog Elastic " << InterpolationPolicy::name() << " Distribution";

    return iss.str();
  }
  static std::string concreteName(
            const Utility::AnalogElasticDistribution<InterpolationPolicy>& instance )
  {
    return name();
  }
};

/*! \brief Type name traits partial specialization for the
 * Utility::UnitAwareAnalogElasticDistribution
 *
 * \details The name function will set the type name that must be used in
 * xml files
 */
template<typename InterpolationPolicy, typename U, typename V>
class TypeNameTraits<Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,U,V> >
{
  public:
  static std::string name()
  {
    std::ostringstream iss;
    iss << "Unit-Aware Analog Elastic " << InterpolationPolicy::name()
        << " Distribution ("
        << Utility::UnitTraits<U>::symbol() << ","
        << Utility::UnitTraits<V>::symbol() << ")";

    return iss.str();
  }
  static std::string concreteName( const Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,U,V>& instance )
  {
    return name();
  }
};

} // end Teuchos namespace

//---------------------------------------------------------------------------//
// Template inludes.
//---------------------------------------------------------------------------//

#include "Utility_AnalogElasticDistribution_def.hpp"

//---------------------------------------------------------------------------//

#endif // end UTILITY_ANALOG_ELASTIC_DISTRIBUTION_HPP

//---------------------------------------------------------------------------//
// end Utility_AnalogElasticDistribution.hpp
//---------------------------------------------------------------------------//
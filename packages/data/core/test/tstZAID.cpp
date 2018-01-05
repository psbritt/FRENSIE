//---------------------------------------------------------------------------//
//!
//! \file   tstZAID.cpp
//! \author Alex Robinson
//! \brief  ZAID class unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <sstream>

// FRENSIE Includes
#include "Data_ZAID.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

typedef std::tuple<
  std::tuple<boost::archive::xml_oarchive,boost::archive::xml_iarchive>,
  std::tuple<boost::archive::text_oarchive,boost::archive::text_iarchive>,
  std::tuple<boost::archive::binary_oarchive,boost::archive::binary_iarchive>,
  std::tuple<Utility::HDF5OArchive,Utility::HDF5IArchive>,
  std::tuple<boost::archive::polymorphic_oarchive*,boost::archive::polymorphic_iarchive*>
  > TestArchives;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a ZAID can be constructed from a raw zaid
FRENSIE_UNIT_TEST( ZAID, raw_zaid_constructor_stable )
{
  Data::ZAID h1_zaid( 1001 );

  FRENSIE_CHECK_EQUAL( h1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h1_zaid.toRaw(), 1001 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_zaid, 1001 );

  Data::ZAID h2_zaid( 1002 );

  FRENSIE_CHECK_EQUAL( h2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h2_zaid.toRaw(), 1002 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_zaid, 1002 );

  Data::ZAID u238_zaid( 92238 );

  FRENSIE_CHECK_EQUAL( u238_zaid.atom(), Data::U_ATOM );
  FRENSIE_CHECK_EQUAL( u238_zaid.atomicNumber(), 92 );
  FRENSIE_CHECK_EQUAL( u238_zaid.atomicMassNumber(), 238 );
  FRENSIE_CHECK_EQUAL( u238_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( u238_zaid.toRaw(), 92238 );
  FRENSIE_CHECK_EQUAL( (unsigned)u238_zaid, 92238 );
}

//---------------------------------------------------------------------------//
// Check that a ZAID can be constructed from a raw zaid
FRENSIE_UNIT_TEST( ZAID, raw_zaid_constructor_metastable )
{
  Data::ZAID h1_meta1_zaid( 1401 );

  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.toRaw(), 1401 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_meta1_zaid, 1401 );

  Data::ZAID h1_meta2_zaid( 1801 );

  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.isomerNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.toRaw(), 1801 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_meta2_zaid, 1801 );

  Data::ZAID h2_meta1_zaid( 1402 );

  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.toRaw(), 1402 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_meta1_zaid, 1402 );

  Data::ZAID h2_meta2_zaid( 1802 );

  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.isomerNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.toRaw(), 1802 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_meta2_zaid, 1802 );

  Data::ZAID u238_meta_zaid( 92638 );

  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atom(), Data::U_ATOM );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atomicNumber(), 92 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atomicMassNumber(), 238 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.toRaw(), 92638 );
  FRENSIE_CHECK_EQUAL( (unsigned)u238_meta_zaid, 92638 );
}

//---------------------------------------------------------------------------//
// Check that a ZAID can be constructed from is components
FRENSIE_UNIT_TEST( ZAID, component_constructor )
{
  Data::ZAID h1_zaid( 1, 1, 0 );

  FRENSIE_CHECK_EQUAL( h1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h1_zaid.toRaw(), 1001 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_zaid, 1001 );

  Data::ZAID h1_meta1_zaid( 1, 1, 1 );

  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid.toRaw(), 1401 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_meta1_zaid, 1401 );

  Data::ZAID h1_meta2_zaid( 1, 1, 2 );

  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.isomerNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid.toRaw(), 1801 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_meta2_zaid, 1801 );

  Data::ZAID h2_zaid( 1, 2, 0 );

  FRENSIE_CHECK_EQUAL( h2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h2_zaid.toRaw(), 1002 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_zaid, 1002 );

  Data::ZAID h2_meta1_zaid( 1, 2, 1 );

  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid.toRaw(), 1402 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_meta1_zaid, 1402 );

  Data::ZAID h2_meta2_zaid( 1, 2, 2 );

  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.atomicMassNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.isomerNumber(), 2 );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid.toRaw(), 1802 );
  FRENSIE_CHECK_EQUAL( (unsigned)h2_meta2_zaid, 1802 );

  Data::ZAID u238_zaid( 92, 238, 0 );

  FRENSIE_CHECK_EQUAL( u238_zaid.atom(), Data::U_ATOM );
  FRENSIE_CHECK_EQUAL( u238_zaid.atomicNumber(), 92 );
  FRENSIE_CHECK_EQUAL( u238_zaid.atomicMassNumber(), 238 );
  FRENSIE_CHECK_EQUAL( u238_zaid.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( u238_zaid.toRaw(), 92238 );
  FRENSIE_CHECK_EQUAL( (unsigned)u238_zaid, 92238 );

  Data::ZAID u238_meta_zaid( 92, 238, 1 );

  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atom(), Data::U_ATOM );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atomicNumber(), 92 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.atomicMassNumber(), 238 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.isomerNumber(), 1 );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid.toRaw(), 92638 );
  FRENSIE_CHECK_EQUAL( (unsigned)u238_meta_zaid, 92638 );
}

//---------------------------------------------------------------------------//
// Check that a ZAID can be copied
FRENSIE_UNIT_TEST( ZAID, copy_constructor )
{
  Data::ZAID h1_zaid( 1, 1, 0 );

  Data::ZAID h1_zaid_copy( h1_zaid );

  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.toRaw(), 1001 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_zaid_copy, 1001 );
  FRENSIE_CHECK_EQUAL( h1_zaid, h1_zaid_copy );
}

//---------------------------------------------------------------------------//
// Check that a ZAID can be assigned
FRENSIE_UNIT_TEST( ZAID, assignment_operator )
{
  Data::ZAID h1_zaid( 1, 1, 0 );

  Data::ZAID h1_zaid_copy = h1_zaid;

  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.atomicMassNumber(), 1 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.isomerNumber(), 0 );
  FRENSIE_CHECK_EQUAL( h1_zaid_copy.toRaw(), 1001 );
  FRENSIE_CHECK_EQUAL( (unsigned)h1_zaid_copy, 1001 );
  FRENSIE_CHECK_EQUAL( h1_zaid,  h1_zaid_copy );
}

//---------------------------------------------------------------------------//
// Check tht a ZAID can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( ZAID, archive, TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_zaid" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    Data::ZAID h1_zaid( 1, 1, 0 );
    Data::ZAID h1_meta1_zaid( 1, 1, 1 );
    Data::ZAID h1_meta2_zaid( 1, 1, 2 );

    Data::ZAID h2_zaid( 1, 2, 0 );
    Data::ZAID h2_meta1_zaid( 1, 2, 1 );
    Data::ZAID h2_meta2_zaid( 1, 2, 2 );

    Data::ZAID u238_zaid( 92, 238, 0 );
    Data::ZAID u238_meta_zaid( 92, 238, 1 );

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h1_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h1_meta1_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h1_meta2_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h2_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h2_meta1_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( h2_meta2_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( u238_zaid ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( u238_meta_zaid ) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  Data::ZAID h1_zaid;
  Data::ZAID h1_meta1_zaid;
  Data::ZAID h1_meta2_zaid;
  
  Data::ZAID h2_zaid;
  Data::ZAID h2_meta1_zaid;
  Data::ZAID h2_meta2_zaid;

  Data::ZAID u238_zaid;
  Data::ZAID u238_meta_zaid;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h1_zaid ) );
  FRENSIE_CHECK_EQUAL( h1_zaid, Data::ZAID( 1, 1, 0 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h1_meta1_zaid ) );
  FRENSIE_CHECK_EQUAL( h1_meta1_zaid, Data::ZAID( 1, 1, 1 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h1_meta2_zaid ) );
  FRENSIE_CHECK_EQUAL( h1_meta2_zaid, Data::ZAID( 1, 1, 2 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h2_zaid ) );
  FRENSIE_CHECK_EQUAL( h2_zaid, Data::ZAID( 1, 2, 0 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h2_meta1_zaid ) );
  FRENSIE_CHECK_EQUAL( h2_meta1_zaid, Data::ZAID( 1, 2, 1 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( h2_meta2_zaid ) );
  FRENSIE_CHECK_EQUAL( h2_meta2_zaid, Data::ZAID( 1, 2, 2 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( u238_zaid ) );
  FRENSIE_CHECK_EQUAL( u238_zaid, Data::ZAID( 92, 238, 0 ) );
  
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( u238_meta_zaid ) );
  FRENSIE_CHECK_EQUAL( u238_meta_zaid, Data::ZAID( 92, 238, 1 ) );
}

//---------------------------------------------------------------------------//
// end tstZAID.cpp
//---------------------------------------------------------------------------//

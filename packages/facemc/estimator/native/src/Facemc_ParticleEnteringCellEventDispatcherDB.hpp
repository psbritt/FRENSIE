//---------------------------------------------------------------------------//
//!
//! \file   Facemc_ParticleEnteringCellEventDispatcherDB.hpp
//! \author Alex Robinson
//! \brief  Particle entering cell event dispatcher database class decl.
//!
//---------------------------------------------------------------------------//

#ifndef FACEMC_PARTICLE_ENTERING_CELL_EVENT_DISPATCHER_DB_HPP
#define FACEMC_PARTICLE_ENTERING_CELL_EVENT_DISPATCHER_DB_HPP

// Boost Includes
#include <boost/unordered_map.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/deref.hpp>

// Teuchos Includes
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "Facemc_ParticleEnteringCellEventDispatcher.hpp"

namespace Facemc{

//! The particle entering cell event dispatcher database class
class ParticleEnteringCellEventDispatcherDB
{

public:

  //! Get the appropriate dispatcher for the given cell id
  static Teuchos::RCP<ParticleEnteringCellEventDispatcher>& getDispatcher(
		    const Geometry::ModuleTraits::InternalCellHandle cell_id );

  //! Attach an observer to the appropriate dispatcher
  static void attachObserver(
		   const Geometry::ModuleTraits::InternalCellHandle cell_id,
		   const ModuleTraits::InternalEstimatorHandle estimator_id,
		   Teuchos::RCP<ParticleEnteringCellEventObserver>& observer );
  
  //! Detach an observer from the appropriate dispatcher
  static void detachObserver(
		    const Geometry::ModuleTraits::InternalCellHandle cell_id,
		    const ModuleTraits::InternalEstimatorHandle estimator_id );

  //! Detach the observer from all dispatchers
  static void detachObserver(
		    const ModuleTraits::InternalEstimatorHandle estimator_id );

  //! Detach all observers
  static void detachAllObservers();

  //! Dispatch the particle entering cell event to the observers
  static void dispatchParticleEnteringCellEvent( 
	      const ParticleState& particle,
	      const Geometry::ModuleTraits::InternalCellHandle cell_entering );

private:

  // Constructor
  ParticleEnteringCellEventDispatcherDB();

  // Typedef for the dispatcher map
  typedef boost::unordered_map<Geometry::ModuleTraits::InternalCellHandle,
			    Teuchos::RCP<ParticleEnteringCellEventDispatcher> >
  DispatcherMap;

  static DispatcherMap master_map;
};

// Dispatch the particle entering cell event to the observers
inline void
ParticleEnteringCellEventDispatcherDB::dispatchParticleEnteringCellEvent( 
	       const ParticleState& particle,
	       const Geometry::ModuleTraits::InternalCellHandle cell_entering )
{
  DispatcherMap::iterator it = 
    ParticleEnteringCellEventDispatcherDB::master_map.find( cell_entering );

  if( it != ParticleEnteringCellEventDispatcherDB::master_map.end() )
    it->second->dispatchParticleEnteringCellEvent( particle, cell_entering );
}

//! Register an observer with the appropriate dispatchers
template<typename Observer, typename EntityHandle>
inline void registerObserver( Teuchos::RCP<Observer>& observer,
			      const Teuchos::Array<EntityHandle>& entity_ids,
			      ParticleEnteringCellEventObserver::EventTag )
{
  // Make sure the Observer class has the corrent event tag
  testStaticPrecondition((boost::is_same<typename boost::mpl::deref<typename boost::mpl::find<typename Observer::EventTags,ParticleEnteringCellEventObserver::EventTag>::type>::type,ParticleEnteringCellEventObserver::EventTag>::type));
  
  Teuchos::RCP<ParticleEnteringCellEventObserver> observer_base = 
    Teuchos::rcp_dynamic_cast<ParticleEnteringCellEventObserver>( observer );
  
  for( unsigned i = 0u; i < entity_ids.size(); ++i )
  {
    ParticleEnteringCellEventDispatcherDB::attachObserver( entity_ids[i],
							   observer->getId(),
							   observer_base );
  }
}

} // end Facemc namespace

#endif // end FACEMC_PARTICLE_ENTERING_CELL_EVENT_DISPATCHER_DP_HPP

//---------------------------------------------------------------------------//
// end Facemc_ParticleEnteringCellEventDispatcherDB.hpp
//---------------------------------------------------------------------------//

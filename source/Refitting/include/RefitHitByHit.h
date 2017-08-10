#ifndef RefitHitByHit_h
#define RefitHitByHit_h 1

#include <marlin/Processor.h>

#include <UTIL/BitField64.h>

#include <EVENT/Track.h>

namespace MarlinTrk{
  class IMarlinTrkSystem;
  class IMarlinTrack;
}


class RefitHitByHit : public marlin::Processor {


public:


  virtual marlin::Processor*  newProcessor() { return new RefitHitByHit ; }

  RefitHitByHit() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( lcio::LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( lcio::LCEvent * evt ) ;


  virtual void check( lcio::LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

protected:

  int FitInit2( Track* track , MarlinTrk::IMarlinTrack* _marlinTrk ) ;

  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection( lcio::LCEvent * evt, std::string colName ) ;

  /* /\* helper function to get relations using try catch block *\/ */
  /* lcio::LCRelationNavigator* GetRelations(lcio::LCEvent * evt, std::string RelName ) ; */

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name ;

  /** output collection name for the not used hits.
   */
  std::string _output_not_used_col_name ;

  /** output track collection name.
   */
  std::string _output_track_col_name ;

  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name ;

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;

  int _n_run ;
  int _n_evt ;

  bool _MSOn=true;
  bool _ElossOn=true;
  bool _SmoothOn=false;
  double _Max_Chi2_Incr=1000.0;
  int _nInitialHits=6;

  float _bField=0.0;

  bool _extrapolateForward=true;

  std::shared_ptr<UTIL::BitField64> _encoder{};


} ;





#endif

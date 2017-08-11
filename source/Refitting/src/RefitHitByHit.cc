#include "RefitHitByHit.h"

#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>

#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/MarlinTrkUtils.h>

#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/BitField64.h>
#include <UTIL/Operators.h>

#include <algorithm>

using namespace lcio;
using namespace marlin;


RefitHitByHit aRefitHitByHit ;


RefitHitByHit::RefitHitByHit() : Processor("RefitHitByHit") {

  // modify processor description
  _description = "RefitHitByHit refits an input  track by creating a small track of the first few hits and then adds individual hits step by step" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACK,
                           "InputTrackCollectionName" ,
                           "Name of the input track collection"  ,
                           _input_track_col_name ,
                           std::string("TruthTracks") ) ;

  registerOutputCollection( LCIO::TRACK,
                            "OutputTrackCollectionName" ,
                            "Name of the output track collection"  ,
                            _output_track_col_name ,
                            std::string("RefittedTracks") ) ;

  registerProcessorParameter("MultipleScatteringOn",
                             "Use MultipleScattering in Fit",
                             _MSOn,
                             bool(true));

  registerProcessorParameter("EnergyLossOn",
                             "Use Energy Loss in Fit",
                             _ElossOn,
                             bool(true));

  registerProcessorParameter("SmoothOn",
                             "Smooth All Mesurement Sites in Fit",
                             _SmoothOn,
                             bool(false));

  registerProcessorParameter("Max_Chi2_Incr",
                             "maximum allowable chi2 increment when moving from one site to another",
                             _Max_Chi2_Incr,
                             _Max_Chi2_Incr);

  registerProcessorParameter("NumberOfHits",
                             "Number of Hits for the initial fit",
                             _nInitialHits,
                             _nInitialHits);

  registerProcessorParameter("extrapolateForward",
                             "if true extrapolation in the forward direction (in-out), otherwise backward (out-in)",
                             _extrapolateForward,
                             _extrapolateForward);

}


void RefitHitByHit::init() {


  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  // usually a good idea to
  printParameters() ;

  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , nullptr , "" ) ;

  ///////////////////////////////

  _encoder = std::make_shared<UTIL::BitField64>( lcio::LCTrackerCellID::encoding_string() );


  if( _trksystem == 0 ){

    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest" )  ) ;

  }

  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;

  _n_run = 0 ;
  _n_evt = 0 ;

}


void RefitHitByHit::processRunHeader( LCRunHeader* ) {

  ++_n_run ;
}

void RefitHitByHit::processEvent( LCEvent * evt ) {

  ++_n_evt ;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection( evt, _input_track_col_name ) ;
  if( not input_track_col ){
    return;
  }

  // establish the track collection that will be created
  LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;
  _encoder->reset();


  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackVec->setFlag( trkFlag.getFlag()  ) ;

  int nTracks = input_track_col->getNumberOfElements()  ;

  streamlog_out(DEBUG4) << " ######### NO OF TRACKS $$$$$$$$$$ " << nTracks << std::endl;

  // loop over the input tracks and refit
  for(int iTrack=0; iTrack< nTracks ; ++iTrack) {

    Track* track = dynamic_cast<Track*>( input_track_col->getElementAt( iTrack ) ) ;

    MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

    EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;

    streamlog_out(DEBUG2) <<"---- tracks n = "<< iTrack << "  n hits = "<<  trkHits.size() << std::endl;

    EVENT::TrackerHitVec::iterator it = trkHits.begin();

    const int nHitsTrack = trkHits.size();
    int nHitsAdded = 0;

    for (int iHit = 0; iHit < _nInitialHits && iHit < nHitsTrack ;++iHit) {
      auto* theHit = trkHits[iHit];
      _encoder->setValue( theHit->getCellID0() ) ;
      const int subdet = (*_encoder)[lcio::LCTrackerCellID::subdet()].value();
      if( subdet > 2 ) break; // only take hits in the vertex barrel or endcap
      marlin_trk->addHit( theHit );
      nHitsAdded++;
    }

    if(nHitsAdded < 3){
      streamlog_out( WARNING ) << "Not enough hits in the vertex endcap/barrel, adding all first hits" << std::endl;
      delete marlin_trk;
      marlin_trk = _trksystem->createTrack();
      for (int iHit = 0; iHit < _nInitialHits && iHit < nHitsTrack ;++iHit) {
        marlin_trk->addHit( trkHits[iHit] );
      }
      nHitsAdded++;
    }

    int init_status = FitInit2(track, marlin_trk) ;

    if (init_status!=0) {
      delete marlin_trk;
      continue;
    }



    streamlog_out(DEBUG4) << "track initialised " << std::endl ;

    int fit_status = marlin_trk->fit();

    if ( fit_status != 0 ){
      delete marlin_trk;
      continue;
    }

    for (int iHit = _nInitialHits; iHit < nHitsTrack; ++iHit) {
      auto nextHit = trkHits[iHit];
      double chi2_increment = 0.;
      bool isSuccessfulFit =
        marlin_trk->addAndFit( nextHit, chi2_increment, _Max_Chi2_Incr ) == MarlinTrk::IMarlinTrack::success;

      streamlog_out(DEBUG4) << " --- hit = " << iHit << std::endl;
      streamlog_out(DEBUG4) << " --- _Max_Chi2_Incr = " << _Max_Chi2_Incr << std::endl;
      streamlog_out(DEBUG4) << " --- increment in the chi2 = " << chi2_increment
                            << "  , max chi2 to accept the hit " << _Max_Chi2_Incr*(iHit+1)  << std::endl;
      streamlog_out(DEBUG4) << " --- isSuccessfulFit = "<< isSuccessfulFit << std::endl ;


    }




    //==============================================================================================================

    IMPL::TrackImpl* lcio_trk = new IMPL::TrackImpl();

    MarlinTrk::IMarlinTrack* marlinTrk = 0 ;

    if( true ) {

      marlinTrk = marlin_trk ;

      bool fit_direction = MarlinTrk::IMarlinTrack::forward ;
      int return_code =  finaliseLCIOTrack( marlin_trk, lcio_trk, trkHits,  fit_direction ) ;

      streamlog_out( DEBUG ) << " *** created finalized LCIO track - return code " << return_code  << std::endl
                             << *lcio_trk << std::endl ;


    }

    // fit finished - get hits in the fit
    std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit* , double> > outliers;

    // remember the hits are ordered in the order in which they were fitted

    marlinTrk->getHitsInFit(hits_in_fit);

    if( hits_in_fit.size() < 3 ) {
      streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  " << trkHits.size() << std::endl;
      delete marlinTrk ;
      delete lcio_trk;
      continue ;
    }

    marlinTrk->getOutliers(outliers);

    std::vector<TrackerHit*> all_hits;
    all_hits.reserve( hits_in_fit.size() + outliers.size() );

    for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }

    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }


    UTIL::BitField64 encoder2(  lcio::LCTrackerCellID::encoding_string() );
    encoder2.reset() ;  // reset to 0
    MarlinTrk::addHitNumbersToTrack(lcio_trk, all_hits, false, encoder2);
    MarlinTrk::addHitNumbersToTrack(lcio_trk, hits_in_fit, true, encoder2);


    streamlog_out( DEBUG4 )  << "RefitHitByHit::processEvent - Hit numbers for track " << lcio_trk->id() << ":  " << std::endl;
    int detID = 0;
    for (size_t ip=0; ip<lcio_trk->subdetectorHitNumbers().size(); ip=ip+2){
      detID++;
      streamlog_out( DEBUG4 )  << "  det id " << detID
                               << " , nhits in track = " << lcio_trk->subdetectorHitNumbers()[ip]
                               << " , nhits in fit = " << lcio_trk->subdetectorHitNumbers()[ip+1]
                               << std::endl;
      if (lcio_trk->subdetectorHitNumbers()[ip] > 0) lcio_trk->setTypeBit( detID ) ;
    }

    trackVec->addElement( lcio_trk );

  } // for loop to the tracks





    //-------------------------------------------------------------------------------------------------------


  evt->addCollection( trackVec , _output_track_col_name ) ;

}


void RefitHitByHit::check( LCEvent* ) {}


void RefitHitByHit::end() {}


LCCollection* RefitHitByHit::GetCollection( LCEvent * evt, std::string colName ){

  LCCollection* col = NULL;

  try{
    col = evt->getCollection( colName.c_str() ) ;
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;
  }

  return col;

}


int RefitHitByHit::FitInit2( Track* track, MarlinTrk::IMarlinTrack* _marlinTrk ){

  TrackStateImpl trackState( TrackState::AtOther,
                             track->getD0(),
                             track->getPhi(),
                             track->getOmega(),
                             track->getZ0(),
                             track->getTanLambda(),
                             track->getCovMatrix(),
                             track->getReferencePoint()
                             );

  bool direction =  _extrapolateForward ? MarlinTrk::IMarlinTrack::forward: MarlinTrk::IMarlinTrack::backward;
  _marlinTrk->initialise( trackState, _bField, direction );

  return MarlinTrk::IMarlinTrack::success ;

}

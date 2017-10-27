#ifndef QuickTrackAssociatorByHitsImpl_h
#define QuickTrackAssociatorByHitsImpl_h

#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "FWCore/Utilities/interface/IndexSet.h"

// Forward declarations
class TrackerHitAssociator;

namespace edm {
 class EDProductGetter;
}

class QuickTrackAssociatorByHitsImpl : public reco::TrackToTrackingParticleAssociatorBaseImpl
{
 public:
  enum SimToRecoDenomType {denomnone,denomsim,denomreco};

  QuickTrackAssociatorByHitsImpl(edm::EDProductGetter const& productGetter,
                                 std::unique_ptr<const TrackerHitAssociator> hitAssoc,
                                 const ClusterTPAssociation *clusterToTPMap,
                                 bool absoluteNumberOfHits,
                                 double qualitySimToReco,
                                 double puritySimToReco,
                                 double cutRecoToSim,
                                 double pixelHitWeight,
                                 bool threeHitTracksAreSpecial,
                                 SimToRecoDenomType simToRecoDenominator);

    reco::RecoToSimCollection associateRecoToSim( const edm::Handle<edm::View<reco::Track> >& trackCollectionHandle,
                                                  const edm::Handle<TrackingParticleCollection>& trackingParticleCollectionHandle) const override;
    reco::SimToRecoCollection associateSimToReco( const edm::Handle<edm::View<reco::Track> >& trackCollectionHandle,
                                                  const edm::Handle<TrackingParticleCollection>& trackingParticleCollectionHandle) const override;
    reco::RecoToSimCollection associateRecoToSim( const edm::RefToBaseVector<reco::Track>& trackCollection,
                                                  const edm::RefVector<TrackingParticleCollection>& trackingParticleCollection) const override;

    reco::SimToRecoCollection associateSimToReco( const edm::RefToBaseVector<reco::Track>& trackCollection,
                                                  const edm::RefVector<TrackingParticleCollection>& trackingParticleCollection) const override;

  //seed
    reco::RecoToSimCollectionSeed associateRecoToSim(const edm::Handle<edm::View<TrajectorySeed> >&,
                                                     const edm::Handle<TrackingParticleCollection>&) const override;

    reco::SimToRecoCollectionSeed associateSimToReco(const edm::Handle<edm::View<TrajectorySeed> >&,
                                                     const edm::Handle<TrackingParticleCollection>&) const override;


 private:
  typedef std::pair<uint32_t,EncodedEventId> SimTrackIdentifiers;

  typedef edm::IndexSet TrackingParticleRefKeySet;

  // - added by S. Sarkar
  static bool tpIntPairGreater(std::pair<edm::Ref<TrackingParticleCollection>,size_t> i, std::pair<edm::Ref<TrackingParticleCollection>,size_t> j) { return (i.first.key()>j.first.key()); }

  template<class T_TrackCollection, class T_TrackingParticleCollection, class T_hitOrClusterAssociator>
  reco::RecoToSimCollection associateRecoToSimImplementation( const T_TrackCollection& trackCollection, const T_TrackingParticleCollection& trackingParticleCollection, const TrackingParticleRefKeySet *trackingParticleKeys, T_hitOrClusterAssociator hitOrClusterAssociator ) const;

  template<class T_TrackCollection, class T_TrackingParticleCollection, class T_hitOrClusterAssociator>
  reco::SimToRecoCollection associateSimToRecoImplementation( const T_TrackCollection& trackCollection, const T_TrackingParticleCollection& trackingParticleCollection, const TrackingParticleRefKeySet *trackingParticleKeys, T_hitOrClusterAssociator hitOrClusterAssociator ) const;


  template<typename T_TPCollection,typename iter> std::vector< std::pair<edm::Ref<TrackingParticleCollection>,double> > associateTrack( const TrackerHitAssociator& hitAssociator, const T_TPCollection& trackingParticles, const TrackingParticleRefKeySet *trackingParticleKeys, iter begin, iter end ) const;
  template<typename T_TPCollection,typename iter> std::vector< std::pair<edm::Ref<TrackingParticleCollection>,double> > associateTrack( const ClusterTPAssociation& clusterToTPMap, const T_TPCollection& trackingParticles, const TrackingParticleRefKeySet *trackingParticleKeys, iter begin, iter end ) const;


  bool trackingParticleContainsIdentifier( const TrackingParticle* pTrackingParticle, const SimTrackIdentifiers& identifier ) const;

  template<typename iter> double getDoubleCount( const TrackerHitAssociator& hitAssociator, iter begin, iter end, TrackingParticleRef associatedTrackingParticle ) const;
  template<typename iter> double getDoubleCount( const ClusterTPAssociation& clusterToTPList, iter begin, iter end, TrackingParticleRef associatedTrackingParticle ) const;

  template<typename iter> std::vector< std::pair<SimTrackIdentifiers,double> > getAllSimTrackIdentifiers( const TrackerHitAssociator& hitAssociator, iter begin, iter end ) const;

  // Added by S. Sarkar
  template<typename iter> std::vector< OmniClusterRef> getMatchedClusters( iter begin, iter end ) const;

  const TrackingRecHit* getHitFromIter(trackingRecHit_iterator iter) const {
    return &(**iter);
  }

  const TrackingRecHit* getHitFromIter(TrackingRecHitCollection::const_iterator iter) const {
    return &(*iter);
  }

  // The last parameter is used to decide whether we cound hits or clusters
  double weightedNumberOfTrackClusters(const reco::Track& track, const TrackerHitAssociator&) const;
  double weightedNumberOfTrackClusters(const TrajectorySeed& seed, const TrackerHitAssociator&) const;
  double weightedNumberOfTrackClusters(const reco::Track& track, const ClusterTPAssociation&) const;
  double weightedNumberOfTrackClusters(const TrajectorySeed& seed, const ClusterTPAssociation&) const;

  // called only by weightedNumberOfTrackClusters(..., ClusterTPAssociation)
  template<typename iter> double weightedNumberOfTrackClusters(iter begin, iter end) const ;

  //void prepareEitherHitAssociatorOrClusterToTPMap( const edm::Event* pEvent, std::unique_ptr<ClusterTPAssociation>& pClusterToTPMap, std::unique_ptr<TrackerHitAssociator>& pHitAssociator ) const;

  edm::EDProductGetter const* productGetter_;
  std::unique_ptr<const TrackerHitAssociator> hitAssociator_;
  const ClusterTPAssociation *clusterToTPMap_;

  double qualitySimToReco_;
  double puritySimToReco_;
  double pixelHitWeight_;
  double cutRecoToSim_;
  SimToRecoDenomType simToRecoDenominator_;
  bool threeHitTracksAreSpecial_;
  bool absoluteNumberOfHits_;

  // Added by S. Sarkar
  //bool useClusterTPAssociation_;
}; // end of the QuickTrackAssociatorByHitsImpl class

#endif // end of ifndef QuickTrackAssociatorByHitsImpl_


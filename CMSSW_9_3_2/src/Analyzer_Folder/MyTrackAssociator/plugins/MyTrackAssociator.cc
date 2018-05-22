// -*- C++ -*-
//
// Package:    SimTracker/TrackAssociatorProducers
// Class:      QuickTrackAssociatorByHitsProducer
//
/**\class QuickTrackAssociatorByHitsProducer QuickTrackAssociatorByHitsProducer.cc SimTracker/TrackAssociatorProducers/plugins/QuickTrackAssociatorByHitsProducer.cc

  Description: Associator to extract data for multilayer perceptron training.

  Implementation:
      For input data please refer to Prof. Dr. Wolfgang Adam
  */
//
// Original Author:  Christoph Bernkopf
//         Created:  Mon, 31 Jul 2017 10:50:34 GMT
//
//

//// ---- Marker 1 ---- ////
// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

//--------------------------------------------------------
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
// if locally saved
//#include "Analyzer_Folder/MyTrackAssociator/interface/QuickTrackAssociatorByHitsImpl.h"
// if used from CMSSW lib
#include "SimTracker/TrackAssociatorProducers/plugins/QuickTrackAssociatorByHitsImpl.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

// Import missing track data
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

// TTree include
#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GsfTrackToVtx imports
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//// ---- Marker 2 ---- ////
//
// class declaration
//
namespace
{
}

class MyTrackAssociator : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit MyTrackAssociator(const edm::ParameterSet &);
    ~MyTrackAssociator();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob() override;
    //    virtual void analyze(edm::StreamID, edm::Event&, const edm::EventSetup&) override;
    virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
    virtual void endJob() override;
    edm::ParameterSet makeHitAssociatorParameters(const edm::ParameterSet &);

    // ----------associator member data ---------------------------
    TrackerHitAssociator::Config trackerHitAssociatorConfig_;
    edm::EDGetTokenT<ClusterTPAssociation> cluster2TPToken_;
    double qualitySimToReco_;
    double puritySimToReco_;
    double pixelHitWeight_;
    double cutRecoToSim_;
    QuickTrackAssociatorByHitsImpl::SimToRecoDenomType simToRecoDenominator_;
    bool threeHitTracksAreSpecial_;
    bool useClusterTPAssociation_;
    bool absoluteNumberOfHits_;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<TrajectorySeed>> TrajectorySeedToken_;
    edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
    edm::EDGetTokenT<edm::View<reco::GsfTrack>> GsfTrackCollectionToken_;
    edm::EDGetTokenT<edm::View<reco::Track>> TrackCollectionToken_;

    // GsfTrackToVtx imports
    edm::ESHandle<GlobalTrackingGeometry> trackingGeometryHandle_;
    edm::ESHandle<MagneticField> magneticFieldHandle_;

    // ----------counting variables ---------------------------
    int indexEvent;
    int assocseedfound;
    int mixturefound;
    double seedsuccessrate;
    double mixturesuccessrate;
    int track_varib_nr;
    int tp_varib_nr;

    // ----------TTree Varibs ---------------------------
    TTree *track_tree;
    float gsf_track[9];
    float seed_assoc_track[9];
    float track_assoc_track[9];
    float tp_track[5];
    float size_nc_weight[3];
    int ic_para[1];
    float stats[5];

    AlgebraicVector5 localPars_;
    AlgebraicMatrix55 localCov_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

MyTrackAssociator::MyTrackAssociator(const edm::ParameterSet &iConfig) : trackerHitAssociatorConfig_(makeHitAssociatorParameters(iConfig), consumesCollector()),
                                                                         qualitySimToReco_(iConfig.getParameter<double>("Quality_SimToReco")),
                                                                         puritySimToReco_(iConfig.getParameter<double>("Purity_SimToReco")),
                                                                         pixelHitWeight_(iConfig.getParameter<double>("PixelHitWeight")),
                                                                         cutRecoToSim_(iConfig.getParameter<double>("Cut_RecoToSim")),
                                                                         threeHitTracksAreSpecial_(iConfig.getParameter<bool>("ThreeHitTracksAreSpecial")),
                                                                         useClusterTPAssociation_(iConfig.getParameter<bool>("useClusterTPAssociation")),
                                                                         absoluteNumberOfHits_(iConfig.getParameter<bool>("AbsoluteNumberOfHits"))
{

    //// ---- Marker 3 ---- ////
    indexEvent = 0;
    assocseedfound = 0;
    mixturefound = 0;
    seedsuccessrate = 0;
    mixturesuccessrate = 0;
    track_varib_nr = 9;
    tp_varib_nr = 5;

    TrajectorySeedToken_ = consumes<edm::View<TrajectorySeed>>(edm::InputTag("electronMergedSeeds"));
    tpToken_ = consumes<TrackingParticleCollection>(edm::InputTag("tpSelection"));
    GsfTrackCollectionToken_ = consumes<edm::View<reco::GsfTrack>>(edm::InputTag("electronGsfTracks"));
    TrackCollectionToken_ = consumes<edm::View<reco::Track>>(edm::InputTag("electronGsfTracks"));

    usesResource("TFileService");
}

MyTrackAssociator::~MyTrackAssociator()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// Set up the parameter set for the hit associator
edm::ParameterSet
MyTrackAssociator::makeHitAssociatorParameters(const edm::ParameterSet &iConfig)
{
    edm::ParameterSet hitAssociatorParameters;
    hitAssociatorParameters.addParameter<bool>("associatePixel", iConfig.getParameter<bool>("associatePixel"));
    hitAssociatorParameters.addParameter<bool>("associateStrip", iConfig.getParameter<bool>("associateStrip"));
    // This is the important one, it stops the hit associator searching through the list of sim hits.
    // I only want to use the hit associator methods that work on the hit IDs (i.e. the uint32_t trackId
    // and the EncodedEventId eventId) so I'm not interested in matching that to the PSimHit objects.
    hitAssociatorParameters.addParameter<bool>("associateRecoTracks", true);
    // add these new ones to allow redirection of inputs:
    hitAssociatorParameters.addParameter<edm::InputTag>("pixelSimLinkSrc", iConfig.getParameter<edm::InputTag>("pixelSimLinkSrc"));
    hitAssociatorParameters.addParameter<edm::InputTag>("stripSimLinkSrc", iConfig.getParameter<edm::InputTag>("stripSimLinkSrc"));

    return hitAssociatorParameters;
}

//// ---- Marker 4 ---- ////
// ------------ method called for each event  ------------
void MyTrackAssociator::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
//MyTrackAssociator::analyze(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    using namespace edm;

    iSetup.get<GlobalTrackingGeometryRecord>().get(trackingGeometryHandle_);
    iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle_);
    MultiTrajectoryStateTransform mtst(&*trackingGeometryHandle_, &*magneticFieldHandle_);

    edm::Handle<edm::View<reco::GsfTrack>> gsfTrackHandle;
    iEvent.getByToken(GsfTrackCollectionToken_, gsfTrackHandle);

    const ClusterTPAssociation *clusterAssoc = nullptr;
    std::unique_ptr<TrackerHitAssociator> trackAssoc;
    trackAssoc = std::make_unique<TrackerHitAssociator>(iEvent, trackerHitAssociatorConfig_);

    edm::Handle<edm::View<TrajectorySeed>> TrajectorySeedHandle;
    iEvent.getByToken(TrajectorySeedToken_, TrajectorySeedHandle);

    edm::Handle<TrackingParticleCollection> tpHandle;
    iEvent.getByToken(tpToken_, tpHandle);

    edm::Handle<edm::View<reco::GsfTrack>> GsfTrackCollectionHandle;
    iEvent.getByToken(GsfTrackCollectionToken_, GsfTrackCollectionHandle);

    edm::Handle<edm::View<reco::Track>> TrackCollectionHandle;
    iEvent.getByToken(TrackCollectionToken_, TrackCollectionHandle);

    std::cout << "\n"
              << "------------------------------------------"
              << "\n"
              << "\n"
              << "--- Output Prints of MyTrackAssociator ---"
              << "\n"
              << "\n"
              << "#GsfTrackCollection = " << GsfTrackCollectionHandle->size() << "\n"
              << "#TrajectorySeeds = " << TrajectorySeedHandle->size() << "\n"
              << "#RecoTracks = " << TrackCollectionHandle->size() << "\n"
              << "#TrackingParticles = " << tpHandle->size() << "\n"
              << std::endl;

    // Associator Function
    auto impl = std::make_unique<QuickTrackAssociatorByHitsImpl>(iEvent.productGetter(),
                                                                 std::move(trackAssoc),
                                                                 clusterAssoc,
                                                                 absoluteNumberOfHits_,
                                                                 qualitySimToReco_,
                                                                 puritySimToReco_,
                                                                 pixelHitWeight_,
                                                                 cutRecoToSim_,
                                                                 threeHitTracksAreSpecial_,
                                                                 simToRecoDenominator_);

    reco::RecoToSimCollectionSeed mySeedToSim = impl->associateRecoToSim(TrajectorySeedHandle, tpHandle);

    reco::RecoToSimCollection myTrackToSim = impl->associateRecoToSim(TrackCollectionHandle, tpHandle);

    //// ---- Marker 5 ---- ////
    // Associator loop
    for (size_t j = 0; j < GsfTrackCollectionHandle->size(); ++j)
    {
        std::cout << "Assocation Loop Number = " << j << std::endl;
        const reco::GsfTrack &gsfTrack = GsfTrackCollectionHandle->at(j);

        ++indexEvent;

        const edm::RefToBase<TrajectorySeed> &mySeedRef = gsfTrack.seedRef();
        reco::RecoToSimCollectionSeed::const_iterator iassocseed = mySeedToSim.find(mySeedRef);

        edm::RefToBase<reco::Track> seedRef(TrackCollectionHandle, j);
        reco::RecoToSimCollection::const_iterator iassoctrack = myTrackToSim.find(seedRef);

        for (int k = 0; k < track_varib_nr; k++)
        {
            gsf_track[k] = 0;
            seed_assoc_track[k] = 0;
            track_assoc_track[k] = 0;
        }

        size_nc_weight[0] = 0;
        size_nc_weight[1] = 0;
        size_nc_weight[2] = 0;
        ic_para[0] = 0;

        for (int k = 0; k < tp_varib_nr; k++)
        {
            tp_track[k] = 0;
        }

        if (iassocseed != mySeedToSim.end())
        {
            size_t kmax = 0;
            double qmax = -1.;

            for (size_t i = 0; i < (*iassocseed).val.size(); i++)
            {

                if ((*iassocseed).val[i].second > qmax)
                {
                    kmax = i;
                    qmax = (*iassocseed).val[i].second;
                }
                else
                {
                }
            }

            if (qmax > 0.)
            {
                const edm::Ref<TrackingParticleCollection> tref_seed = (*iassocseed).val[kmax].first;
                ++assocseedfound;
            }
            else
            {
                std::cout << "No sim to reco seed!"
                          << "\n"
                          << std::endl;
            }
        }

        //// ---- Marker 6 ---- ////
        if (iassoctrack != myTrackToSim.end())
        {

            std::cout << "Sim to reco track found!" << std::endl;
            std::cout << "#myTrackToSim Size = " << myTrackToSim.size() << std::endl;
            std::cout << "(*iassoctrack).val.size() = " << (*iassoctrack).val.size() << std::endl;

            size_t kmax = 0;
            double qmax = -1.;

            for (size_t i = 0; i < (*iassoctrack).val.size(); i++)
            {
                std::cout << "Get correct seed Loop = " << i << std::endl;

                if ((*iassoctrack).val[i].second > qmax)
                {
                    kmax = i;
                    qmax = (*iassoctrack).val[i].second;
                }
                else
                {
                    std::cout << "iassoctrack not filled, bad quality!" << std::endl;
                }
            }

            if (qmax > 0.)
            {
                std::cout << "qmax < 0 found!" << std::endl;

                const edm::Ref<TrackingParticleCollection> tref_track = (*iassoctrack).val[kmax].first;
                std::cout << "track_assoc_track writen!" << std::endl;

                const reco::GsfTrack *igsf = &(*gsfTrackHandle)[j];

                TrajectoryStateOnSurface innTSOS = mtst.innerStateOnSurface(*igsf);
                TrajectoryStateOnSurface vtxTSOS = mtst.extrapolatedState(innTSOS, GlobalPoint(tref_track->vx(), tref_track->vy(), tref_track->vz()));

                LocalVector assocp(vtxTSOS.surface().toLocal(GlobalVector(tref_track->px(), tref_track->py(), tref_track->pz())));
                LocalPoint assocv(vtxTSOS.surface().toLocal(GlobalPoint(tref_track->vx(), tref_track->vy(), tref_track->vz())));

                if ((assocv.x() * assocv.x() + assocv.y() * assocv.y()) < 0.001)
                {

                    size_nc_weight[0] = j;
                    size_nc_weight[1] = vtxTSOS.components().size();
                    size_nc_weight[2] = 1;
                    ic_para[0] = -1;

                    localPars_ = vtxTSOS.localParameters().vector();
                    localCov_ = vtxTSOS.localError().matrix();

                    tp_track[0] = tref_track->charge() / tref_track->p();
                    tp_track[1] = assocp.x() / assocp.z();
                    tp_track[2] = assocp.y() / assocp.z();
                    tp_track[3] = assocv.x();
                    tp_track[4] = assocv.y();

                    //// ---- Marker 7 ---- ////
                    track_tree->Fill();

                    std::cout << "Gsf total mixture filled" << std::endl;

                    for (size_t ic = 0; ic < vtxTSOS.components().size(); ++ic)
                    {
                        //// ---- Marker 8 ---- ////
                        LocalVector assocp(vtxTSOS.surface().toLocal(GlobalVector(tref_track->px(), tref_track->py(), tref_track->pz())));
                        LocalPoint assocv(vtxTSOS.surface().toLocal(GlobalPoint(tref_track->vx(), tref_track->vy(), tref_track->vz())));

                        //// ---- Marker 9 ---- ////
                        size_nc_weight[0] = j;
                        size_nc_weight[1] = vtxTSOS.components().size();
                        size_nc_weight[2] = vtxTSOS.components()[ic].weight();
                        ic_para[0] = ic;

                        //                              std::cout << "size = " << typeid(vtxTSOS.components().size()).name() << "ic " << typeid(ic).name() << std::endl;

                        localPars_ = vtxTSOS.components()[ic].localParameters().vector();
                        localCov_ = vtxTSOS.components()[ic].localError().matrix();

                        tp_track[0] = tref_track->charge() / tref_track->p();
                        tp_track[1] = assocp.x() / assocp.z();
                        tp_track[2] = -assocp.y() / assocp.z();
                        tp_track[3] = assocv.x();
                        tp_track[4] = assocv.y();

                        track_tree->Fill();

                        std::cout << "Loop " << ic << " complete & filled" << std::endl;
                    }

                    ++mixturefound;
                    std::cout << "mixturefound # increased!"
                              << "\n"
                }
                else
                {
                    std::cout << "Gsf Vertex too far from origin, no Gaussian Mix produced"
                              << "\n"
                              << std::endl;
                }
            }
        }
        else
        {
            std::cout << "No sim to reco track!"
                      << "\n"
                      << std::endl;
        }
    }
    //// ---- Marker 10 ---- ///

    seedsuccessrate = float(assocseedfound) / float(indexEvent);
    mixturesuccessrate = float(mixturefound) / float(indexEvent);

    std::cout << "indexEvent = " << indexEvent << "\n"
              << "mixturefound = " << mixturefound << "\n"
              << "p found = " << mixturesuccessrate << "\n"
              << "\n"
              << "------------------------------------------"
              << "\n"
              << std::endl;

    for (int k = 0; k < 5; k++)
    {
        stats[k] = 0;
    }

    stats[0] = indexEvent;
    stats[3] = mixturefound;
    stats[4] = mixturesuccessrate;
}

//------------------------------------------------------------------------------

// ------------ method called once each job just before starting event loop  ------------
void MyTrackAssociator::beginJob()
{

    using namespace edm;

    // initialize tree
    edm::Service<TFileService> fs;
    track_tree = fs->make<TTree>("track_associator_tree", "Associator tree with branches");
    track_tree->Branch("size_nc_weight", &size_nc_weight, "size_nc_weight[3]/F");
    track_tree->Branch("ic_para", &ic_para, "ic_para[1]/I");
    track_tree->Branch("localPars", &localPars_);
    track_tree->Branch("localCov", &localCov_);
    track_tree->Branch("tp_track", &tp_track, "tp_track[5]/F");
    track_tree->Branch("stats", &stats, "stats[5]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void MyTrackAssociator::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyTrackAssociator::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTrackAssociator);

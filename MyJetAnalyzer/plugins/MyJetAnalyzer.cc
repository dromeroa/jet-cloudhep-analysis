//*****************************************************************************************
//*****************************************************************************************
//// -*- C++ -*-
//
// Package:    AnalizaJet/MyJetAnalyzer
// Class:      MyJetAnalyzer
// 
/**\class MyJetAnalyzer MyJetAnalyzer.cc AnalizaJet/MyJetAnalyzer/plugins/MyJetAnalyzer.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]*/
//
// Original Author: David Romero Abad  
//         Created:  Wed, 16 Apr 2025 19:05:17 GMT
//***************************************************************************************
//***************************************************************************************


//****************************************************************************
//*********************         LIBRERIAS        *****************************
//****************************************************************************
//Librerias del sistema
#include <memory>
#include <TMath.h>
// Librerias del usuario
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"
#include "TH1.h"
//Clases extraida de la informacion de jets
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//clases para guardar los datos 
#include "TTree.h"
#include "TFile.h"
#include<vector>
#include "TRandom3.h"


//****************************************************************************
//*********************    CLASE MyJetAnalzer    *****************************
//****************************************************************************
// Para recordar el concepto de clase en C++, se puede revisar: https://www.w3schools.com/cpp/cpp_classes.asp
class MyJetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MyJetAnalyzer(const edm::ParameterSet&);
      virtual ~MyJetAnalyzer();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
     // virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // declarar la etiqueta (token) de entrada para la coleccion de jets, rho y vertices primarios. 
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::VertexCollection> pvToken_;

      // ----------member data ---------------------------
      // Variables para la correcion de la energia de jets (JEC)
      //JEC (Jet Energy Corrections) son factores aplicados para corregir la energía medida de los jets
      //Estas correciones necesitan de archivos externos de entrada con parametros calibrados: Ejemplo: Summer20UL18_V5_MC_L2Relative_AK4PFchs.txt
      std::vector<std::string> jecPayloadNames_;
      std::string              jecL2_; // corrige variaciones en la respuesta del detector como función de la pseudorapidez (η).(Para Datos y MC)
      std::string              jecL3_; // Compensa efectos como la no linealidad del detector y pérdidas de energía no instrumentadas (pT)(Para Datos y MC)
      std::string              jecRes_;// correcion residual:(Datos) Compensa las diferencias residuales, basándose en comparaciones de datos con simulaciones
      std::string              jetJECUncName_;
      std::string              jetResName_;
      std::string              sfName_;
      boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
      boost::shared_ptr<FactorizedJetCorrector> jecL2L3_;
      bool isData;

      // Variables para las correciones por Jet Energy resolution (JER)
      // JER es una medida de cuán precisa es la reconstrucción del momento transversal (pT) del jet comparado con el valor real a nivel de partícula
      // Las correciones tambien necesitan de archivos externos, ejemplo: Summer20UL18_JRV3_MC_PtResolution_AK4PFchs.txt
      JME::JetResolution resolution;
      JME::JetResolutionScaleFactor resolution_sf;

      // Definimos el Tree  https://root.cern/manual/trees/
      // Definimos los vectores donde vamos a guardar las colecciones
      // Variables para el fatjet
      int numfatjet; 
      TTree *mtree;
      std::vector<float> fatjet_e;
      std::vector<float> fatjet_pt;
      std::vector<float> fatjet_eta;
      std::vector<float> fatjet_phi;
      std::vector<float> fatjet_ch;
      std::vector<float> fatjet_mass;
      std::vector<float> fatjet_corrpt;
      std::vector<float> fatjet_corrptUp;
      std::vector<float> fatjet_corrptDown;
      std::vector<float> fatjet_corrptSmearUp;
      std::vector<float> fatjet_corrptSmearDown;
      std::vector<float> fatjet_corrmass;
      std::vector<float> fatjet_corre;
      std::vector<float> fatjet_corrpx;
      std::vector<float> fatjet_corrpy;
      std::vector<float> fatjet_corrpz;
      std::vector<float> fatjet_prunedmass;
      std::vector<float> fatjet_softdropmass;
      std::vector<float> fatjet_tau1;
      std::vector<float> fatjet_tau2;
      std::vector<float> fatjet_tau3;
      std::vector<float> fatjet_subjet1btag;
      std::vector<float> fatjet_subjet2btag;
      std::vector<float> fatjet_subjet1hflav;
      std::vector<float> fatjet_subjet2hflav;
};

//****************************************************************************************************************************
//***********************************************     INICIALIZAR LAS CLASES        ******************************************
//****************************************************************************************************************************
//Declara e inicializa varias variables miembro privadas de la clase con valores que vienen del archivo de configuración Python.

MyJetAnalyzer::MyJetAnalyzer(const edm::ParameterSet& iConfig):
  // token del jet que se conecta con el archivo de configuracion python
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  // token de los PV
  pvToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  // token de las correciones del jet aplicadas a Monte Carlo
  jecL2_(iConfig.getParameter<edm::FileInPath>("jecL2Name").fullPath()), // JEC level payloads for test            
  jecL3_(iConfig.getParameter<edm::FileInPath>("jecL3Name").fullPath()), // JEC level payloads for test            
  jecRes_(iConfig.getParameter<edm::FileInPath>("jecResName").fullPath()),
  jetJECUncName_(iConfig.getParameter<edm::FileInPath>("jetJECUncName").fullPath()), // JEC uncertainties
  jetResName_(iConfig.getParameter<edm::FileInPath>("jerResName").fullPath()), // JER Resolutions
  sfName_(iConfig.getParameter<edm::FileInPath>("jerSFName").fullPath()), // JER Resolutions
  isData(iConfig.getParameter<bool>("isData"))
{
  //crea una instancia de servicio para manejar la creación y almacenamiento de un archivo ROOT.
  edm::Service<TFileService> fs;
  //crea un árbol TTree dentro de un archivo ROOT utilizando el servicio fs
  mtree = fs->make<TTree>("JetsTree", "JetsTree");

   // Prepara y carga las JEC necesarias para aplicar correcciones a los jets, dependiendo de si estás procesando datos reales o simulación (MC). 
   jecPayloadNames_.push_back(jecL2_);
   jecPayloadNames_.push_back(jecL3_);
   if( isData == true ) jecPayloadNames_.push_back(jecRes_);
   std::vector<JetCorrectorParameters> vParL2L3;
   for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNames_.begin(),
          payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
          JetCorrectorParameters pars(*ipayload);
          vParL2L3.push_back(pars);
   }
                                                                                            
  // Este código configura el módulo para: Corregir jets usando factores calibrados (L2+L3).Evaluar incertidumbres en esa corrección.Acceder a información sobre la resolución de jets. Leer variables del evento necesarias como rho y vértices primarios.

  jecL2L3_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vParL2L3) );
  jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jetJECUncName_) );

  resolution = JME::JetResolution(jetResName_);
  resolution_sf = JME::JetResolutionScaleFactor(sfName_);

  edm::InputTag rhotag("fixedGridRhoFastjetAll");
  rhoToken_ = consumes<double>(rhotag);
  //edm::InputTag pvtag("offlineSlimmedPrimaryVertices");
  //pvToken_ = consumes<reco::VertexCollection>(pvtag);

// Se definen las ramas del Tree y sus nombres
  mtree->Branch("numberfatjet",&numfatjet);
  mtree->GetBranch("numberfatjet")->SetTitle("Number of Fatjets");
  mtree->Branch("fatjet_e",&fatjet_e);
  mtree->GetBranch("fatjet_e")->SetTitle("Uncorrected Fatjet energy");
  mtree->Branch("fatjet_pt",&fatjet_pt);
  mtree->GetBranch("fatjet_pt")->SetTitle("Uncorrected Transverse Fatjet Momentum");
  mtree->Branch("fatjet_eta",&fatjet_eta);
  mtree->GetBranch("fatjet_eta")->SetTitle("Fatjet Eta");
  mtree->Branch("fatjet_phi",&fatjet_phi);
  mtree->GetBranch("fatjet_phi")->SetTitle("Fatjet Phi");
  mtree->Branch("fatjet_ch",&fatjet_ch);
  mtree->GetBranch("fatjet_ch")->SetTitle("Fatjet Charge");
  mtree->Branch("fatjet_mass",&fatjet_mass);
  mtree->GetBranch("fatjet_mass")->SetTitle("Fatjet Mass");
  mtree->Branch("fatjet_prunedmass",&fatjet_prunedmass);
  mtree->GetBranch("fatjet_prunedmass")->SetTitle("L2+L3-corrected pruned mass of Fatjet");
  mtree->Branch("fatjet_softdropmass",&fatjet_softdropmass);
  mtree->GetBranch("fatjet_softdropmass")->SetTitle("L2+L3-corrected softdrop mass of Fatjet");
  mtree->Branch("fatjet_tau1",&fatjet_tau1);
  mtree->GetBranch("fatjet_tau1")->SetTitle("N-subjettiness tau_1 of Fatjet");
  mtree->Branch("fatjet_tau2",&fatjet_tau2);
  mtree->GetBranch("fatjet_tau2")->SetTitle("N-subjettiness tau_2 of Fatjet");
  mtree->Branch("fatjet_tau3",&fatjet_tau3);
  mtree->GetBranch("fatjet_tau3")->SetTitle("N-subjettiness tau_3 of Fatjet");
  mtree->Branch("fatjet_subjet1btag",&fatjet_subjet1btag);
  mtree->GetBranch("fatjet_subjet1btag")->SetTitle("Leading softdrop subjet 1 b discriminant");
  mtree->Branch("fatjet_subjet2btag",&fatjet_subjet2btag);
  mtree->Branch("fatjet_subjet1hflav",&fatjet_subjet1hflav);
  mtree->GetBranch("fatjet_subjet1hflav")->SetTitle("Leading softdrop subjet 1 hadron flavour");
  mtree->Branch("fatjet_subjet2hflav",&fatjet_subjet2hflav);
  mtree->GetBranch("fatjet_subjet2hflav")->SetTitle("Leading softdrop subjet 2 hadron flavour");
}

MyJetAnalyzer::~MyJetAnalyzer()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//****************************************************************************************************************************
//***********************************************  ANALYZER - PARTE PRINCIPAL DEL PROGRAMA    ********************************
//****************************************************************************************************************************
// ------------ method called for each event  ------------
void
MyJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 //Para evitar escribir edm y std cada vez. Sin esto tendriamos que escribir: edm::Handle<pat::JetCollection> fatjets; y std::vector<float> fatjet_pt;	
  using namespace edm;
  using namespace std;

  // Este bloque de código extrae del evento tres insumos clave para un análisis de jets:La colección de fat jets (por ejemplo, AK8 jets para estudios de bosones o top). El valor de rho, que representa la densidad de energía en el evento. A menudo se utiliza para la corrección del pileup en la calibración de la energía de los jets o en el aislamiento usado para correcciones por pileup. Los vértices primarios, usados para identificar el origen de los objetos físicos en el evento. Para esto se utilizan tokens previamente definidos: fatjetToken_, rhoToken_, pvToken_.

  //fatjets
  Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);

  // Densidad de energia (rho)
  Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);

  // vertices
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(pvToken_, vertices); 

  // vertices
  //Handle<reco::VertexCollection> vertices;
  //iEvent.getByToken(pvToken_, vertices);

  // Este bloque asegura que los vectores estén vacíos antes de llenarlos con la información de jets del nuevo evento. Esto es esencial para evitar resultados incorrectos en el análisis.

  numfatjet = 0;
  fatjet_e.clear();
  fatjet_pt.clear();
  fatjet_eta.clear();
  fatjet_phi.clear();
  fatjet_ch.clear();
  fatjet_mass.clear();
  fatjet_corrpt.clear();
  fatjet_corrptUp.clear();
  fatjet_corrptDown.clear();
  fatjet_corrptSmearUp.clear();
  fatjet_corrptSmearDown.clear();
  fatjet_corrmass.clear();
  fatjet_corre.clear();
  fatjet_corrpx.clear();
  fatjet_corrpy.clear();
  fatjet_corrpz.clear();
  fatjet_prunedmass.clear();
  fatjet_softdropmass.clear();
  fatjet_tau1.clear();
  fatjet_tau2.clear();
  fatjet_tau3.clear();
  fatjet_subjet1btag.clear();
  fatjet_subjet2btag.clear();
  fatjet_subjet1hflav.clear();
  fatjet_subjet2hflav.clear();

  // Definimos las variables de correcion
  double corrpt;
  double corrUp, corrDown;
  float ptscale, ptscale_down, ptscale_up;

  // Fijamos el valor minimo del pT en 200 GeV
  int min_pt = 200;
   
  // verifica si el Handle llamado fatjets contiene una colección válida obtenida del evento
  // Antes de usar la colección, es buena práctica (y a veces necesario) verificar si el Handle es válido,
  if(fatjets.isValid()){

      //Recorre cada fat jet en la colección fatjets.
      for (const pat::Jet &fatjet : *fatjets){

      // Obtiene la versión no corregida (o con corrección de nivel 0) del jet actual.	      
      pat::Jet uncorrFatjet = fatjet.correctedJet(0);

      //Asigna el valor del momento transversal corregido del jet (pt) a la variable corrpt. 
      corrpt = fatjet.pt();

      // Inicializa los factores de corrección hacia arriba (corrUp) y hacia abajo (corrDown) en 1.0 (sin corrección por el momento).
      corrUp = 1.0;
      corrDown = 1.0;

      // Si el valor absoluto de (eta) es menor que 5, se usa su valor directamente. Si no, se limita a 4.99 (para evitar problemas numéricos)
      if( fabs(fatjet.eta()) < 5) jecUnc_->setJetEta( fatjet.eta() );
      else jecUnc_->setJetEta( 4.99 );

      // Se establece el valor del momento transversal corregido del jet como entrada para calcular la incertidumbre.
      jecUnc_->setJetPt( corrpt );
 
      // Se calcula el factor de corrección hacia arriba: 1 más la magnitud de la incertidumbre (modo ascendente, con argumento 1).
      corrUp = (1 + fabs(jecUnc_->getUncertainty(1)));

      // Mismo procedimiento para el caso descendente
      if( fabs(fatjet.eta()) < 5) jecUnc_->setJetEta( fatjet.eta() );
      else jecUnc_->setJetEta( 4.99 );
      jecUnc_->setJetPt( corrpt );
      corrDown = (1 - fabs(jecUnc_->getUncertainty(-1)));

      ptscale = 1;
      ptscale_down = 1;
      ptscale_up = 1;
 
      // Si es simulacion Monte Carlo (MC) 
      if(!isData) {

	//Se define un conjunto de parámetros para calcular la resolución JER. corrpt es el pt corregido del jet. fatjet.eta() es la pseudorrapidez del jet.
	//rho es la densidad de energía del evento (usada en varias correcciones).

        JME::JetParameters JERparameters = {{JME::Binning::JetPt, corrpt}, {JME::Binning::JetEta, fatjet.eta()}, {JME::Binning::Rho, *(rhoHandle.product())}};
 
        //Calcula la resolución JER nominal y los scale factors (centrales y variaciones hacia arriba y abajo). Estos valores se usan para hacer el smearing del pt del jet. 

	float res = resolution.getResolution(JERparameters);
        float sf = resolution_sf.getScaleFactor(JERparameters);
        float sf_up = resolution_sf.getScaleFactor(JERparameters, Variation::UP);
        float sf_down = resolution_sf.getScaleFactor(JERparameters, Variation::DOWN);

	// fatjet es un jet reconstruido en el evento. fatjet.genJet() obtiene el jet generado (GenJet) correspondiente a este fatjet. Los jets generados son los jets que están presentes en la simulación del proceso físico (MC), es decir, son los jets producidos en la simulación del modelo de física que se está estudiando. El puntero genFatjet será nullptr si no se encuentra un jet generado correspondiente.

	// deltaR < 0.4: Esto significa que los jets generados y reconstruidos deben estar cercanos en el espacio de 4-vectores, es decir, que deben estar bien emparejados en términos de su dirección. deltaPt <= 3 * corrpt * res: Esta condición asegura que la diferencia en pt no sea demasiado grande. Si el jet generado y el jet corregido tienen una gran discrepancia en pt, no se realiza el smearing.

        const reco::GenJet *genFatjet = fatjet.genJet();
        bool smeared = false;
        if(genFatjet){
          double deltaPt = fabs(genFatjet->pt() - corrpt);
          double deltaR = reco::deltaR(genFatjet->p4(),fatjet.p4());
          if ((deltaR < 0.4) && deltaPt <= 3*corrpt*res){
            ptscale = max(0.0, 1 + (sf - 1.0)*(corrpt - genFatjet->pt())/corrpt);
            ptscale_down = max(0.0, 1 + (sf_down - 1.0)*(corrpt - genFatjet->pt())/corrpt);
            ptscale_up = max(0.0, 1 + (sf_up - 1.0)*(corrpt - genFatjet->pt())/corrpt);
            smeared = true;
          }
        }

	// Este fragmento de código está relacionado con el proceso de smearing de jets en simulaciones de física de partículas. En este caso, se aplica un ajuste al pt de un fatjet (un tipo de jet con un pt elevado), teniendo en cuenta la resolución de energía del jet (JER) y los factores de escala.

	// Este código tiene como objetivo aplicar un ajuste aleatorio (smearing) a la energía de los jets para simular la incertidumbre y resolución de las mediciones en los experimentos reales. Este proceso se basa en la resolución (res) y los factores de escala (sf, sf_down, sf_up) que representan la variabilidad y las incertidumbres en los datos simulados.


        if (!smeared) {
        TRandom3 JERrand;

        JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
        ptscale = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf*sf - 1.0)));

        JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
        ptscale_down = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_down*sf_down - 1.0)));

        JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
        ptscale_up = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_up*sf_up - 1.0)));
          }
        }

	// Este fragmento de código realiza un análisis detallado de un jet, aplicando un ajuste estocástico (smearing) a sus propiedades de energía, momento y masa. Se almacenan tanto las propiedades originales como las corregidas, junto con las variaciones de las correcciones para estudiar la incertidumbre en el proceso.

        if( ptscale*corrpt <= min_pt) continue;

        pat::Jet smearedFatjet = fatjet;
        smearedFatjet.scaleEnergy(ptscale);

        fatjet_e.push_back(uncorrFatjet.energy());
        fatjet_pt.push_back(uncorrFatjet.pt());
        fatjet_eta.push_back(uncorrFatjet.eta());
        fatjet_phi.push_back(uncorrFatjet.phi());
        fatjet_ch.push_back(uncorrFatjet.charge());
        fatjet_mass.push_back(uncorrFatjet.mass());

        fatjet_corrpt.push_back(smearedFatjet.pt());
        fatjet_corrptUp.push_back(corrUp*smearedFatjet.pt());
        fatjet_corrptDown.push_back(corrDown*smearedFatjet.pt());
        fatjet_corrptSmearUp.push_back(ptscale_up*smearedFatjet.pt()/ptscale);
        fatjet_corrptSmearDown.push_back(ptscale_down*smearedFatjet.pt()/ptscale);
        fatjet_corrmass.push_back(smearedFatjet.mass());
        fatjet_corre.push_back(smearedFatjet.energy());
        fatjet_corrpx.push_back(smearedFatjet.px());
        fatjet_corrpy.push_back(smearedFatjet.py());
        fatjet_corrpz.push_back(smearedFatjet.pz());
 
	double corrL2L3 = 1;
        jecL2L3_->setJetEta( uncorrFatjet.eta() );
        jecL2L3_->setJetPt ( uncorrFatjet.pt() );
        jecL2L3_->setJetE  ( uncorrFatjet.energy() );
        jecL2L3_->setJetA  ( smearedFatjet.jetArea() );
        jecL2L3_->setRho   ( *(rhoHandle.product()) );
        jecL2L3_->setNPV   ( vertices->size() );
        corrL2L3 = jecL2L3_->getCorrection();

        fatjet_prunedmass.push_back(corrL2L3*(double)smearedFatjet.userFloat("ak8PFJetsCHSPrunedMass"));
        fatjet_softdropmass.push_back(corrL2L3*(double)smearedFatjet.userFloat("ak8PFJetsCHSSoftDropMass"));
        fatjet_tau1.push_back((double)smearedFatjet.userFloat("NjettinessAK8:tau1"));
        fatjet_tau2.push_back((double)smearedFatjet.userFloat("NjettinessAK8:tau2"));
        fatjet_tau3.push_back((double)smearedFatjet.userFloat("NjettinessAK8:tau3"));

        auto const & sdSubjets = smearedFatjet.subjets("SoftDrop");
        int nSDSubJets = sdSubjets.size();
        if(nSDSubJets > 0){
          pat::Jet subjet1 = sdSubjets.at(0);
          fatjet_subjet1btag.push_back(subjet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
          fatjet_subjet1hflav.push_back(subjet1.hadronFlavour());
        }else{
          fatjet_subjet1btag.push_back(-999);
          fatjet_subjet1hflav.push_back(-999);
        }
        if(nSDSubJets > 1){
          pat::Jet subjet2 = sdSubjets.at(1);
          fatjet_subjet2btag.push_back(subjet2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
          fatjet_subjet2hflav.push_back(subjet2.hadronFlavour());
        }else{
          fatjet_subjet2btag.push_back(-999);
          fatjet_subjet2hflav.push_back(-999);
        }

        ++numfatjet;
    }
  }

  mtree->Fill();
  return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyJetAnalyzer::endJob() 
{
}


// ------------ method called when starting to processes a run  ------------
void
MyJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
MyJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
MyJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MyJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}




// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyJetAnalyzer);

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing ('analysis')
isData = False
if len(sys.argv) > 2:
    try:
        isData = eval(sys.argv[2])
        sys.argv.pop( 2 )
        print "isData is set to ",isData
    except:
        pass
options.parseArguments()
isMC = True
if isData: isMC = False

# Punto de partida en un archivo del configuracion de CMSSW 
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"

# Selecciona el numero maximo de eventos (con -1 corre sobre todoa los eventos)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Configuracion necesaria para manejar tracks transitorios si se requiere, cargar la configuracion de la geometria del detector y cargar la configuracion del campo magnetico del detector CMS 
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Define el archivo de datos de entrada
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:0634456A-08C2-E511-A0C1-001E6739722E.root'
    )
)


if isData:
    process.source.fileNames = cms.untracked.vstring(
    'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleElectron/MINIAOD/08Jun2016-v1/10000/001A703B-B52E-E611-BA13-0025905A60B6.root'
        )

    #---- Aplicar el filtro de calidad de datos usando el archivo JSON. Este ejemplo es para datos de 2015
    #---- Debe hacerse despues de la definicion de process.source
    #---- Asegurate de que la ubicacion del archivo coincida con tu configuracion
    goodJSON = "data/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt"
    myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


#---- Estas dos lineas son necesarias si necesitas acceso a la base de datos de condiciones. Por ejemplo, para obtener correcciones de energia de jets, prescales de triggers, etc.
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#---- Descomenta y adapta una linea como esta si estas accediendo a la base de datos de condiciones mediante archivos de instantaneas (snapshots) de CVMFS 
#     (requiere tener instalado el cliente de CVMFS)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
#---- Si el contenedor ya tiene archivos de base de datos locales disponibles, descomenta lineas como las de abajo 
#     en lugar de las correspondientes lineas anteriores
if isData: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
else: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_mcRun2_asymptotic_RunIIFall15DR76_v1.db')
#---- The global tag must correspond to the needed epoch (comment out if no conditions needed)
if isData: process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
else: process.GlobalTag.globaltag = "76X_mcRun2_asymptotic_RunIIFall15DR76_v1"


process.mypvertex = cms.EDAnalyzer('VertexAnalyzer',
                                   vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   beams=cms.InputTag("offlineBeamSpot"))

process.mygenparticle = cms.EDAnalyzer('GenParticleAnalyzer',
                                       pruned=cms.InputTag("prunedGenParticles"),
                                       #---- Collect particles with specific "status:pdgid"
                                       #---- if 0:0, collect them all 
                                       input_particle = cms.vstring("1:11","1:13","1:22","2:15"))


#---Se comienza con las correciones del jets (JEC)

JecString = 'MC'
if isData: JecString = 'DATA'

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets

#------- Aplicamos el identificador de los filtros de ruido para los fatjet

process.looseAK8Jets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                    filterParams = pfJetIDSelector.clone(),
                                    src = cms.InputTag("slimmedJetsAK8"))


#----------- Aplicamos las correciones de energia de los jets del 2015
process.patJetCorrFactorsReapplyJECAK8 = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("looseAK8Jets"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
        payload = 'AK8PFchs'
        )
if isData: process.patJetCorrFactorsReapplyJECAK8.levels.append('L2L3Residual')
process.slimmedJetsAK8NewJEC = updatedPatJets.clone(
    jetSource = cms.InputTag("looseAK8Jets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8")),
)


#--------- Configura el analizador de fatjets (Aqui se llama a los archivos txt)
process.myfatjets = cms.EDAnalyzer('FatjetAnalyzer',
                fatjets = cms.InputTag("slimmedJetsAK8NewJEC"),
                isData = cms.bool(isData),
                                jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L2Relative_AK8PFchs.txt'),
                                jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L3Absolute_AK8PFchs.txt'),
                                jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt'),
                jetJECUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt'),
                                jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_PtResolution_AK8PFchs.txt'),
                                jerSFName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'), # AK8 == AK4
                )


#----------- Trigger
#-----------Revisar cuales Trigger son mas compatibles con nuestro analisis
process.mytriggers = cms.EDAnalyzer('TriggerAnalyzer',
                              processName = cms.string("HLT"),
                              #---- These are example triggers for 2012
                              #---- Wildcards * and ? are accepted (with usual meanings)
                               #---- If left empty, all triggers will run              
                              triggerPatterns = cms.vstring("HLT_L2DoubleMu23_NoVertex_v*","HLT_Mu12_v*", "HLT_Photon20_CaloIdVL_v*", "HLT_Ele22_CaloIdL_CaloIsoVL_v*", "HLT_Jet370_NoJetID_v*"),
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")
                              )



#----- RUN THE JOB! -----#
process.TFileService = cms.Service("TFileService", fileName=cms.string("Nuevo-FatJets-tree.root"))

if isData:
    process.p = cms.Path(process.mypvertex+ process.looseAK8Jets+process.patJetCorrFactorsReapplyJECAK8+process.slimmedJetsAK8NewJEC+process.myfatjets
                     )
else:
    process.p = cms.Path( process.mypvertex+process.mygenparticle+ process.looseAK8Jets+process.patJetCorrFactorsReapplyJECAK8+process.slimmedJetsAK8NewJEC+process.myfatjets
                     )
process.maxEvents.input = options.maxEvents
process.TFileService.fileName = options.outputFile
if len(options.inputFiles) > 0:
    process.source.fileNames=options.inputFiles

print "Processing for maxEvents =  ",process.maxEvents.input
print "Processing input files "
for fl in process.source.fileNames:
    print "  > ",fl
print "Output filename : ",process.TFileService.fileName                

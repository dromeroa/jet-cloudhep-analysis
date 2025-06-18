# Proyecto Analyze - Jets
## Para clonar el repositorio

```git clone https://github.com/dromeroa/jet-cloudhep-analysis```


```cd MyJetAnalyzer```

## Para correr el código

```cmsenv```

```scram b```

```cd python```

```cmsRun MC_ConfFile_cfg.py outputFile=TreeFatjets.root```


* Usar un archivo root de MC para correr el programa:

```0634456A-08C2-E511-A0C1-001E6739722E.root```

ubicado en el drive


## Para graficar

Una vez obtenido el archivo root: ```TreeFatjets.root```

En la misma carpeta, ejecutar:

```root plot_Nsub.C```


## Tambien se puede usar este dataset para correr (TTbar)

```https://opendata.cern.ch/record/19978```



## Para copiar archivos desde mi espacio en la nube al Docker

```docker cp /home/ubuntu/David/cms_open_data_work/CMSSW_7_6_7/src/jet-cloudhep-analysis/MyJetAnalyzer/python my_od:/code/CMSSW_7_6_7/src/jet-cloudhep-analysis/MyJetAnalyzer/python```

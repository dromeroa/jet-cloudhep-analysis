#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void plot_Nsub() {
    // Abrir el archivo ROOT que contiene el TTree
    TFile *inputFile1 = TFile::Open("TreeFatjets.root", "READ");
    if (!inputFile1 || inputFile1->IsZombie()) {
        std::cerr << "Error abriendo el archivo ROOT fatjet.root ." << std::endl;
        return;
    }

    // Obtener el TTree
    TTree *tree1 = (TTree*)inputFile1->Get("myfatjets/Events");
    if (!tree1) {
        std::cerr << "No se encontró el TTree 'myfatjets/Events' en TreeFatjets.root." << std::endl;
        inputFile1->Close();
        return;
    }

    // Crear punteros para los vectores de cada TTree
    std::vector<float> *tau1 = nullptr;
    std::vector<float> *tau2 = nullptr;
    std::vector<float> *tau3 = nullptr;

    // Establecer las ramas en el TTree para los vectores
    tree1->SetBranchAddress("fatjet_tau1", &tau1); // Solo necesitas pasar los punteros a los vectores
    tree1->SetBranchAddress("fatjet_tau2", &tau2);
    tree1->SetBranchAddress("fatjet_tau3", &tau3);
 

    // Crear un histograma para almacenar los resultados
    TH1F *h1 = new TH1F("histograma", "Histograma de tau21", 100, 0, 1);
    TH1F *h2 = new TH1F("histograma", "Histograma de tau31", 100, 0, 1); 


    // Leer los eventos y llenar el histograma con los elementos de tau1
    Long64_t nEntries = tree1->GetEntries();  // Usar Long64_t para el número de entradas
    std::cout << "Número de entradas en el TTree: " << nEntries << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree1->GetEntry(i);  // Cargar el evento i

        // Verificar que tau1 no sea null antes de usarlo
        if (tau1) {
            for (size_t j = 0; j < tau1->size(); ++j) {
                // Llenar el histograma con los valores de tau1
                h1->Fill(tau2->at(j)/tau1->at(j));
                h2->Fill(tau3->at(j)/tau1->at(j));
            }
        }
    }

    // Crear un canvas para dibujar el histograma
    TCanvas *c1 = new TCanvas("canvas", "Histograma de tau21", 800, 600);
    h1->SetAxisRange(0, 1, "X");
    h1->Draw();
    c1->SaveAs("tau21.pdf");


    TCanvas *c2 = new TCanvas("canvas", "Histograma de tau31", 800, 600);
    h2->SetAxisRange(0, 1, "X");
    h2->Draw();
    c2->SaveAs("tau31.pdf");








    // Cerrar el archivo ROOT
    inputFile1->Close();
}


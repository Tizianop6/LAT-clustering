double fitfunc(double* x, double* p) {
   //funzione di fit
   double fitval;

   double norm = p[0];
   fitval=p[0]; //+ p[1]*x[0];

   return fitval;

}

void fitfolding() {

   const int mpt = 300000;
   double x1[mpt];
   double iniz=2.39562e+08,fin=6.49606e+08, ampiezza=410044624.;
   int npt = 0;
   int npt2=0;
   double folding=0;
   // read data file
   ifstream in;
   in.open("tempicluster"); //file con i MET di arrivo di ogni cluster, quello fornito è un esempio

   while ( kTRUE ) {
      in >> x1[npt];
      if ( ! in.good() ) break;
      npt++;

   }
   in.close();

   cout << npt;
   cout << endl << x1[0] << endl << x1[npt-1];
   in.close();

   printf(" found %d points\n", npt);

   TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 500, 500); //canvas di root

   TH1D* h1 = new TH1D("h1", "Istogramma Dati Raccolti", 544,2.39562e+08, 6.49606e+08); //estremi su cui fare l'istogramma
    h1->SetFillColor(kBlue);
   for ( int i = 0; i < npt; i++ ) {

     h1->Fill(x1[i]);

   }
   
   TF1* f1 = new TF1("f1", fitfunc, -100.00, 100.0, 2); //questa in realtà non serve, non verrà usata dopo
   TF1 g1("g1", "pol0", 0., 10e8); //si dichiara la funzione di fit
   g1.SetParameter(1.,2.); // parametri iniziali (ne chiede due per qualche misterioso motivo)
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptFit(kTRUE);

   c1->cd();
   h1->Fit("g1");
   ofstream myfile;
   myfile.open("50bin.txt", ios::out | ios::trunc);
   for (float i=1; i< 200; i=i+0.1){
   cout << endl<< i<< endl << endl<< ampiezza/(i*3600.*24.) << endl;
   folding=i*3600.*24.;
   
   	   TH1D* h1 = new TH1D("h1", "Istogramma Dati Raccolti", 50,0., folding);
   	      for ( int j = 0; j < npt; j++ ) {

		     h1->Fill(fmod(x1[j],folding)); //si riempie l'istogramma con i dati foldati
		   }

   	   h1->Fit("g1"); //fit 
   	   myfile << i << "	"<< g1.GetChisquare()/g1.GetNDF() << "	" << g1.GetNDF() << "\n"; //si dà in output il tempo di folding, il chi^2 ridotto, i gradi di libertà 
   }
   myfile.close();
   f1->SetLineColor(kRed);
   c1->cd();
   g1.Draw("same");

}

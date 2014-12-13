#include "junctionModeler.h"
#include <QApplication>
#include <QDebug>
#include <iostream>
#include <algorithm>

using namespace std;
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    JunctionModeler w;

    //Send argmument count and argument values to the window ()
    w.argc = argc;
    for(int i = 0; i < argc; i++){
        w.argv.append(argv[i]);
    }
    if(std::find(w.argv.begin(), w.argv.end(), "--tutorial") != w.argv.end()) {
        cout << "Welcome to Junction Modeler Windowless Mode! " << endl;
        cout << "Please read carefully the instructions to follow. " << endl;
        w.tutorial();
    }else if (std::find(w.argv.begin(), w.argv.end(), "--help") != w.argv.end()){
        cout << ":: Getting Started with Junction Modeler :: " << endl;
        cout << "For tutorial please enter the following : " << endl;
        cout << "path_to_modeler.exe --tutorial <'export_location'> <filename> " << endl;
        cout << "Example : " << endl;
        cout << "path_to_modeler.exe --tutorial 'C:/Users/Joe/Desktop' data " << endl;
        cout << "For standard mode enter the following : " << endl;
        cout << "path_to_modeler.exe 'export_location' filename --flag <value> " << endl;
        cout << "Flag Options :  " << endl;
        cout << "--vin : input voltage . Accepts double [in volts]" << endl;
        cout << "--r_auto : auto R to amplify to 2.5 * input voltage, accepts nothing." << endl;
        cout << "--rin : input resistance. Accepts double [in ohms]" << endl;
        cout << "--threshold<1..4> : IV thresholds 1 and 2 are enabled by default. ie --threshold1 0.001 " << endl;
        cout << "--account_for_hysteresis : Enables --thresholds3 and --thresholds4 , accepts nothing." << endl;
        cout << "--freqin : Input frequency : accepts double [in Hz]" << endl;
        cout << "--terms_per_fit : Nterms/IV fit routine, accepts integer" << endl;
        cout << "--max_t : max time axis of waveform (for calculating FFT) accepts double" << endl;
        cout << "--manual_mode : program expects manual A and B coefficients to be inputted instead of IV Curve, accepts nothing" << endl;
        cout << "--manualA<1..2> : A coeff's of IV curve, ex. --manualA2 0.000001 , accepts double" << endl;
        cout << "--manualB<1..2> : B coeff's of IV curve, ex. --manualB1 1 , accepts double" << endl;
        cout << "--exportLocation : path to where file is to be exported, path set to desktop by default, accepts string" << endl;
        cout << "--exportFolderName : folder created to store all the files generated from the application, accepts string" << endl;
        cout << "--tutorial : Goes through guided mode with extra hints" << endl;
        cout << "--guided : Goes through each input individually instead of passing as flags" << endl;
        cout << "--show_defaults : Shows the default values for each flag if not provided" << endl;
    }else if (std::find(w.argv.begin(), w.argv.end(), "--show_defaults") != w.argv.end()){
        for (int i=0;i<18;i++){
            cout << "--" << w.set_arguments.arg_name[i] << " : " << w.set_arguments.arg_vals[i] << endl;
        }
    } else if (argc > 1){
        char icing[] = "-";

        for (unsigned int i = 0; i < 18; ++i)
        {
           w.argv[i].erase (std::remove(w.argv[i].begin(), w.argv[i].end(), icing[0]), w.argv[i].end());
           cout << w.argv[i] << endl;
        }

    }else {
        w.show();
    }

    return a.exec();
}

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QClipboard>
#include <QVector>
#include <QDebug>

namespace Ui {
class MainWindow;
}

class JunctionModeler : public QMainWindow
{
    Q_OBJECT

public:
    // Get argumnets incase ran as Excel Plugin or as console / command based application
    int argc;
    QVector<std::string> argv;


    void tutorial();
    explicit JunctionModeler(QWidget *parent = 0);
    QVector<double> xdata;
    QVector<double> ydata;
    QVector<double> output_x_data;
    QVector<double> output_y_data;
    QVector<double> waveform_x_data;
    QVector<double> waveform_y_data;
    QVector<double> peaks_x_data;
    QVector<double> peaks_y_data;
    QVector<double> power_x_data;
    QVector<double> power_y_data;
    bool accountForHysteresis = false;
    bool inManualMode = false;
    bool useOldData = false;
    ~JunctionModeler();

    template <typename T, size_t N> const T* mybegin(const T (&a)[N]) { return a; }
    template <typename T, size_t N> const T* myend  (const T (&a)[N]) { return a+N; }

    static double absolute(double x){
        if (x < 0)
            return x*-1;
        else
            return x;
    }
    static double sgn(double x){
        double y = (x >= 0) - (x < 0);
        return y;
    }
    struct Hash_Key {
        std::string arg_name[18] = {"vin","rin","r_auto","threshold1","threshold2","threshold3","threshold4","account_for_hysteresis",
                                 "freqin","terms_per_fit","max_t","manual_mode","manualA1","manualA2","manualB1",
                                 "manualB2", "exportLocation","exportFolderName"};
        std::string arg_vals[18] = {"0.3","1000","true","0.001","0.001","0.001","0.001","false",
                                 "150","280","0.226667","false","0.000001","0.000001","1",
                                 "1", "desktop","data"};
       // typedef enum {vin = 0,rin=1,r_auto=2,threshold1=3,threhold2=4,threshold3=5,threshold4=6,account_for_hysteresis=7,freqin=8,
       //              terms_per_fit=9,max_t=10,manual_mode=11,manualA1=12,manualA2=13,manualB1=14,manualB2=15,exportLocation=16,
       //              exportFolderName=17} arg_codes;
    }set_arguments;

    struct VectorSet {
        QVector<double> x;
        QVector<double> y;
        double x_max = 0;
    };
    struct hyperbolic_model {
        QVector<double> x;
        QVector<double> y;
        double A;
        double B;
        double minBound;
        double maxBound;
        bool errors = false;
        double best_error = 0;
    };
    struct SlicedVector {
        VectorSet upwards_positive;
        VectorSet downwards_positive;
        VectorSet downwards_negative;
        VectorSet upwards_negative;
        bool errors = false;
    };

    struct VoltagePointer {
        double input;
        double output;
    };
    struct CustomPoint {
        double min;
        double max;
    };
    struct RootObject {
        double error;
        double value;
        double index;
        QVector<double> parent;
        CustomPoint range;
    };

    static double CubicInterpolate(
       double y0,double y1,
       double y2,double y3,
       double mu)
    {
       double a0,a1,a2,a3,mu2;

       mu2 = mu*mu;
       a0 = y3 - y2 - y0 + y1;
       a1 = y0 - y1 - a0;
       a2 = y2 - y0;
       a3 = y1;

       return(a0*mu*mu2+a1*mu2+a2*mu+a3);
    }

    static int pow2roundup (int x)
    {
        if (x < 0)
            return 0;
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        x |= x >> 32;
        x |= x >> 64;
        return x+1;
    }

    Ui::MainWindow *ui;
private slots:
    void on_pushButton_clicked();

    void on_checkBox_clicked(bool checked);

    void on_checkBox_2_clicked(bool checked);

    void on_rescale_power_clicked();

    void on_rescale_power2_clicked();

    void on_exportIV_clicked();

    void on_exportPeaks_clicked();

    void on_exportWaveform_clicked();

    void on_exportPower_clicked();

    void on_rescale_waveform_clicked();

    void on_rescale_waveform2_clicked();

    void on_exportSettings_clicked();

    void on_exportAll_clicked();

    void on_useOldData_clicked();

//private:


};

#endif // MAINWINDOW_H

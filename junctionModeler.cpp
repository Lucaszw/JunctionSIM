#include "junctionModeler.h"
#include "ui_mainwindow.h"
#include <iostream>
#include "stdio.h"
#include <armadillo>
#include <QClipboard>
#include <QString>
#include <QDebug>
#include <QMovie>
#include <QFile>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_symbol.h>
#include <qwt_legend.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include "persistence1d.hpp"
#include "rootfinder.h"
#include <QErrorMessage>
#include <QStandardPaths>
#include <QDir>

using namespace arma;
using namespace std;


JunctionModeler::JunctionModeler(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    // setup labels for all input params
    // as well as there default values
    ui->setupUi(this);

}

JunctionModeler::~JunctionModeler()
{
    delete ui;
}

void exportToDesktop(QString filename, QVector<double> xdata, QVector<double> ydata, QString prefix, bool inFolder, QString filepath = "desktop") {

    QFile file;
    if (filepath == "desktop")
        filepath = QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);

    if (inFolder){
        if (!QDir(filepath + "\\" +prefix).exists())
            QDir().mkdir(filepath + "\\" +prefix);
        file.setFileName(filepath + "\\" +prefix + "\\"+ prefix + "_" + filename + ".txt");
    }else {
        file.setFileName(filepath + "\\"+ prefix + "_" + filename + ".txt");
    }
    if (!file.open(QIODevice::WriteOnly)) {
        std::cerr << "Cannot open file for writing: "
                  << qPrintable(file.errorString()) << std::endl;
        return;
    }

    QTextStream out(&file);
    for(int i=0;i< xdata.size();i++)
        out << xdata[i] << "\t" << ydata[i]  << endl;
    qDebug() << file.fileName();

    file.close();
}

void sortSet(double* xdata, double* ydata, int n){
    // sorts data without forgeting which x point belongs
    // to which y point.

    for(int i=0; i<n; i++){
        for(int j=0; j<n-1; j++){
            if(xdata[j]>xdata[j+1]){
                double temp_x = xdata[j+1];
                double temp_y = ydata[j+1];
                xdata[j+1] = xdata[j];
                ydata[j+1] = ydata[j];
                xdata[j] = temp_x;
                ydata[j] = temp_y;
            }
        }
    }
    return;

}



JunctionModeler::VectorSet grabFrequencies(QVector<double> x, QVector<double> y, double err = 0){
    // Find peaks of FFT spectrum given by a certain threshold
    QVector<float> xfft;
    QVector<float> yfft;
    for (int i=0;i<x.size();i++){
        xfft.append((float)x[i]);
        yfft.append((float)y[i]);
    }
    // Make sure that your xdata is linearly spaced
    fvec xfft_vec = fvec(xfft.toStdVector());
    vector<float> yfft_vec;
    double xmax = xfft_vec.max();
    double xmin = xfft_vec.min();
    vector<float> xnew = conv_to< vector<float> >::from(linspace(xmin,xmax,xfft_vec.size()));

    // Solve for a y value at every point of xnew (as all are evenly spaced
    // interpolate about the midpoint (mu = 0.5)
    for (int i=3;i<yfft.size();i++)
        yfft_vec.push_back( (float) JunctionModeler::CubicInterpolate(yfft[i-3],yfft[i-2],yfft[i-1],yfft[i],0.5));

    p1d::Persistence1D cubic_extrema;

    cubic_extrema.RunPersistence(yfft.toStdVector());
    vector<int> indx_min;
    vector<int> indx_max;
    cubic_extrema.GetExtremaIndices(indx_min, indx_max,err,false);


    cout << indx_min.size();

    JunctionModeler::VectorSet max_indecies;

    // sort y values from smallest to largest
    for (int i=0;i<indx_min.size();i++)
        max_indecies.y << y[indx_max[i]];
    qSort(max_indecies.y.begin(),max_indecies.y.end());

    // from i to size, look at the index, then loop through y, and find where it is now placed.
    for (int i=0;i<max_indecies.y.size();i++){
        int k;
        for (int j=0;j<max_indecies.y.size();j++){
            if (max_indecies.y[i] == y[indx_max[j]]){
                k = j;
                break;
            }
        }
        // save corresponding x:
        max_indecies.x.append(x[indx_max[k]]);
    }

    // get the maximum x value (used for scaling the x-axis of the power spectra)
    QVector<double> x_temp = max_indecies.x;
    qSort(x_temp.begin(),x_temp.end());

    max_indecies.x_max = x_temp.last();

    return max_indecies;

}

JunctionModeler::CustomPoint range(QVector<double> myvector, bool in_abs = false){
    JunctionModeler::CustomPoint range;
    range.min = myvector[0];
    range.max = myvector[0];

    for (int i=0;i<myvector.size();i++){
        if (in_abs){
            myvector[i] = JunctionModeler::absolute(myvector[i]);
            range.min = JunctionModeler::absolute(range.min);
            range.max = JunctionModeler::absolute(range.max);
        }
        if (range.min > myvector[i])
            range.min = myvector[i];
        if (range.max < myvector[i])
            range.max = myvector[i];
    }
    return range;
}


JunctionModeler::SlicedVector slice_data(JunctionModeler::VectorSet inputData,double threshold,bool accountForHysteresis = false, bool fourSweepMode = true){

    /* Finds general direction and quadrant of data and slices accordingly
        Code by : Lucas Zeer-Wanklyn. Property of NINT */

    // The data is sorted into buckets, the top and bottom of each bucket can be filled with bad data,
    // but the middle should always be the most clean using this algorithm
    // There are 3 primary (and one optional) steps to this function

    // 1. Gather the slope and direction roughly, and put the data in corresponding buckets
    // Once filled, we should expect to not have 100% accuracy however this sorting method
    // will ensure that both : the bucket mostly contains data corresponding to a single set
    // the noise in the data is minimalized close to the mid point of the set
    // Furthermore, the number of items in the set, roughly estimates the number of points
    // of the set (given that all buckets contain roughly the same amount of incorrect data

    // 2. If we look at the index values of all sets, vs. the index values of a single set
    // we know that there is a linear relationship
    // thus linearly space the data as integers between the bounds of the bucket and the midpoint
    // of the buckets set (as this is the most clean)

    // 3. After linearly spacing the data, the items should overlapp slightly on the edges
    // this is due to noisy data at the edges of all sets, and as we can assume that
    // the noise is normaly distributed between all sets, therefore the median of the
    // start and end points of each bucket is roughly the correct starting point for each set

    // 4. If the results are not satisfying, a threshold can be applied as an offset in favor
    // of the one of the neighbouring buckets. I.E. instead of stretching to the median, one bucket
    // can stretch more to add more data points to itself.



    // Find the average step size in the x-direction to avoid dividing by zero, and
    // to help with directional analysis



    int data_range[2] = {(int)inputData.x.size()/2,inputData.x.size()-1};
    JunctionModeler::SlicedVector sliced_data;
    // If the user does not have enough infor for a full IV curve , then may as
    // well sort the data (IE no distinction between up sweep and down sweep)
    if (!fourSweepMode){
          for (int i=0;i<data_range[1] + 1; i++){
              inputData.x[i] *= -1;
              inputData.y[i] *= -1;
          }
        // ensure that the x data is linearly spaced
        vec temp_x_data = linspace(inputData.x[0],inputData.x[data_range[1]],data_range[1]+1);
        double del_x = temp_x_data(1)-temp_x_data(0);
        if (del_x == 0){
            sliced_data.errors = true;
            return sliced_data;
        }
        // get position where data becomes positive:
        double n_switch = -1*temp_x_data(0)/del_x;

        double * xdata[2];
        double * ydata[2];

        xdata[0] = &inputData.x[0];
        ydata[0] = &inputData.y[0];
        xdata[1] = &inputData.x[n_switch];
        ydata[1] = &inputData.y[n_switch];
        sortSet(xdata[0],ydata[0],inputData.x.size()/2 -1);
        sortSet(xdata[1],ydata[1],inputData.x.size()/2 -1);

        accountForHysteresis = false;
    }



    double abs_del_x;

    if (inputData.x.size() < 100 || inputData.y.size() < 100){
        sliced_data.errors = true;
        qDebug() << "Bad data!";
        return sliced_data;
    }
    abs_del_x = JunctionModeler::absolute(JunctionModeler::absolute(inputData.x[3])-JunctionModeler::absolute(inputData.x[0]))/( 3 * 4);

    for (int i=0;i<2;i++){
        abs_del_x += JunctionModeler::absolute(JunctionModeler::absolute(inputData.x[data_range[i]] ) - JunctionModeler::absolute(inputData.x[data_range[i] - 3 ]))/(3 * 4) ;
    }


    QVector<double> buckets[4];
    for (int i=2;i<inputData.y.size();i++){

        // find the direction the curve is moving in and its region accross the y-axis of the curve:
        double slope_i = ( inputData.y[i] - inputData.y[i-2] )/(2*abs_del_x);
        double region_i = inputData.y[i-2] + inputData.y[i-1] + inputData.y[i];

        // normalize :
        slope_i = slope_i/abs(slope_i);
        region_i = region_i/abs(region_i);

        // There are four possible cases for where the point may belong to
        // they can be distinguished based on slope + 2*region after normalization:
        // this formula is arbitrary, but is an algebraic represnetation of combining the two normalized sets:
        // if both sets only contain normalized data, and one set is weighted double then , the order of operations
        // can produce only four possible outcomes +1 + +2, +1 + -2, -1 + -2, and -1 + +2

        switch (int(slope_i + 2*region_i)) {
          case 3:
                buckets[0].append(i);
            break;
          case 1:
                buckets[1].append(i);
            break;
          case -3:
                buckets[2].append(i);
            break;
          case -1:
                buckets[3].append(i);
            break;
          default:
                qDebug() << "Strange data with slope = " << slope_i << "and region = " << 2*region_i;
            break;
        }
    }

    int starting_bounds[4];
    int ending_bounds[4];

    // Get the bounds of the data sets based on extrapolation from the midpoint:

    for (int i=0;i<4;i++){

        if (buckets[i].size() < 15){
            // Impossible to divide data into four regions
            // check if the data can be split into two regions instead of 4
            if ((buckets[0].size() + buckets[2].size()) < 15 || (buckets[1].size() + buckets[3].size()) < 15) {
                sliced_data.errors = true;
                qDebug() << "Bad data!";
                return sliced_data;
            }else {
                buckets[i].append(0);
                buckets[i].append(0);
                buckets[i].append(0);
            }

        }

        int mid = buckets[i].size()/2;
        double b = buckets[i][mid] - mid;
        starting_bounds[i] = b;
        ending_bounds[i] = buckets[i].size() + b;
    }

    // This method is great for removing unwanted data around the edges, but unfortunately can lead to overlapping quadrants
    // to fix this find the mean distance between two buckets, and stretch the edges to match it's bounding means
    int bucket_order[4] = {0,1,2,3};
    int temp;

    for(int i=0;i<4;i++){
        for(int j=i+1;j<4;j++){
            if(ending_bounds[i]>ending_bounds[j]){
                temp=i;
                bucket_order[i]=bucket_order[j];
                bucket_order[j]=temp;
            }
        }
    }

    int new_edges[4][2];
    new_edges[0][0] = 0; new_edges[0][1] = int((ending_bounds[bucket_order[0]]+starting_bounds[bucket_order[1]])/2);
    new_edges[1][0] = new_edges[0][1]+1; new_edges[1][1] = int((ending_bounds[bucket_order[1]]+starting_bounds[bucket_order[2]])/2);
    new_edges[2][0] = new_edges[1][1]+1; new_edges[2][1] = int((ending_bounds[bucket_order[2]]+starting_bounds[bucket_order[3]])/2);
    new_edges[3][0] = new_edges[2][1]+1; new_edges[3][1] = buckets[3][buckets[3].size()-1];


    // Return the data to the user sliced into the four quadrants:
    for (int i=0;i<inputData.x.size();i++){
        if (i <= new_edges[0][1] && i >= new_edges[0][0]){
            sliced_data.upwards_positive.x.append(inputData.x[i]);
            sliced_data.upwards_positive.y.append(inputData.y[i]);
        }
        else if ( i <= new_edges[1][1] && i >= new_edges[1][0]){
            if (accountForHysteresis){
                sliced_data.downwards_positive.x.append(inputData.x[i]);
                sliced_data.downwards_positive.y.append(inputData.y[i]);
            }else {
                sliced_data.upwards_positive.x.append(inputData.x[i]);
                sliced_data.upwards_positive.y.append(inputData.y[i]);
            }
        }
        else if (i <= new_edges[2][1] && i >= new_edges[2][0]) {
            sliced_data.downwards_negative.x.append(inputData.x[i]);
            sliced_data.downwards_negative.y.append(inputData.y[i]);
        }
        else if ( i <= new_edges[3][1] && i >= new_edges[3][0]){
            if (accountForHysteresis){
                sliced_data.upwards_negative.x.append(inputData.x[i]);
                sliced_data.upwards_negative.y.append(inputData.y[i]);
            } else {
                sliced_data.downwards_negative.x.append(inputData.x[i]);
                sliced_data.downwards_negative.y.append(inputData.y[i]);
            }
        }
    }
    if ((sliced_data.upwards_negative.x.empty() && sliced_data.downwards_negative.x.empty())
         || (sliced_data.upwards_positive.x.empty() && sliced_data.downwards_positive.x.empty())){
        sliced_data.errors = true;
        qDebug() << "Bad data!";
        return sliced_data;
    }

    if (!accountForHysteresis){
        double * xdata; double * xdata2;
        double * ydata; double * ydata2;

        xdata = &sliced_data.upwards_positive.x[0];
        ydata = &sliced_data.upwards_positive.y[0];
        xdata2 = &sliced_data.downwards_negative.x[0];
        ydata2 = &sliced_data.downwards_negative.y[0];

        // If not accounting for hysteresis combine the two sets of data together for the both positive and negative x directions
        sortSet(xdata,ydata,sliced_data.upwards_positive.x.size());
        sortSet(xdata2,ydata2,sliced_data.downwards_negative.x.size());

    }

    return sliced_data;

}


JunctionModeler::hyperbolic_model fit_sinh(JunctionModeler::VectorSet inputData, double threshold, bool manual_mode= false, double A_manual = 0, double B_manual = 0, double N_terms = 280){

    // Fit data to A*sinh(b*x)
    // By Lucas Zeer-Wanklyn, Property of NINT

    /* Run fitting algorithm if not in manula mode,
        else skip, and use the coefficients given by user*/

    JunctionModeler::CustomPoint bounds;

    double A_best = pow(1,-3);
    double B_best = 0.1;
    double best_error = 0;

    if (manual_mode == false){

        // Crop the data to optimize center region
        // This is found via the minimum cost given by the variable "threshold"
        // Get the slope of the data for neighbouring points,
        // set the maximum error point to be at the bounds (since these points have
        // the least effect on the junction's current )

        // If the cost to start or end at said point is less than a percentage given
        // of the maximum error, crop the data. This optimization method assumes that
        // the center of the input data is more important than the outer points,
        // one can vary the fitting by change the weights of the importance of the outer
        // points from freshold (ie , increasing threshold => decreasing weights of outer
        // region, decreasing threshold => increasing weights of outer region

        // A value for threshold > 1 means that the outer regions are more important
        // (and should be avoided unless unable to acheive a good fit otherwise)

        // get the spacing between two points in the vector
        // the data may not be evenly spaced thus it is better to find the average distance
        // between two points



        // re-sort data from highest to lowest given the data is sliced

        if (JunctionModeler::absolute(inputData.x[inputData.x.size()-1]) <= JunctionModeler::absolute(inputData.x[0])){
            double temp_set_x[inputData.x.size()];
            double temp_set_y[inputData.y.size()];
            // Get dx, if negative reverse data:

            for (int i=0;i<inputData.x.size();i++){
                temp_set_x[i] = inputData.x[inputData.x.size()-i-1];
                temp_set_y[i] = inputData.y[inputData.x.size()-i-1];
            }

            for (int i=0;i<inputData.x.size();i++){
                inputData.x[i] =  temp_set_x[i];
                inputData.y[i] =  temp_set_y[i];
            }
        }


        double dx = 0;
        for (int i=0;i<inputData.x.size();i++)
            dx += inputData.x[i];
        dx /= inputData.x.size();

        // defince a struct to hold derivative data for the data manipulation method:
        struct deriv_holder {
            // require three points, one centered and two edges
            double x_values[3];
            double y_values[3];
            double dx;
            double getslope(){
                double myslope;
                for(int i=1;i<3;i++){
                    double yv = (y_values[i] - y_values[i-1])/(2*(dx));
                    if ( yv != yv)
                        yv = 1000;
                    myslope += yv;
                }
                return myslope;
            }
        }dx_data;
        // after initializing the deriv struct, set the data spacing:
        dx_data.dx = dx;

        double maxerror = 0;

        // get the maximum expected error, by taking the derivative of the
        // data at the bounds.
        // if one were to optimize based on one point alone,
        // you should expect the maximum difference in your fit to be
        // here:

        for (int i=1;i<3;i++){
            for (int j=i-1; j<i+2; j++){
                dx_data.x_values[j] = inputData.x[j];
                dx_data.y_values[j] = inputData.y[j];
            }
            maxerror -= JunctionModeler::absolute(dx_data.getslope());
            maxerror *= -1;
        }

        // The weighted positions are dependent on the bounds for which
        // the error is below the given threshold
        bounds.min = 0;
        bounds.max = inputData.x.size()-1;
        // loop through array strating 3 from the beginning to ensure accurate slope values
        // try to look for two bounds (ie look for when the slope decreases below a certain
        // threshold and vice vera

        // as we are unaware of the slope range it is better to start assuming a small error threshold
        for (int i=4;i<inputData.x.size();i++){
            // get three points (one to left and one to right, and one centered at j)
            for (int j=0; j<2; j++){
                dx_data.x_values[j] = inputData.x[i - 2 + j];
                dx_data.y_values[j] = inputData.y[i - 2 + j];
            }
            // see how much your data is changing :
            double oldslope = dx_data.getslope();
            double newslope = dx_data.getslope();
            double error = JunctionModeler::absolute(newslope) - JunctionModeler::absolute(oldslope);

            // check if your data is different enough from the edges (ie
            // is X * more important

            if ( error < maxerror*threshold)
                bounds.min = i;
            else if (error >= maxerror*threshold && bounds.min > 0){
                bounds.max = i;
                break;
            }
        }

        // Fit data between 4 points usually gives an accurate approximation
        // step forward, steb backward,
        // step forward again (but to a smaller degree) then backward...
        // then increase steps again...

        double pivot_point[14];
        // Could probably put path into a loop,
        // but is easy two write as points as isn't too complicated :
        pivot_point[0] = bounds.max;
        pivot_point[1] = bounds.min;
        pivot_point[2] = floor(bounds.min + (bounds.max - bounds.min)/1.2);
        pivot_point[3] = floor(bounds.min + (bounds.max - bounds.min)/1.5);
        pivot_point[4] = floor(bounds.min + (bounds.max - bounds.min)/2);
        pivot_point[5] = floor(bounds.min + (bounds.max - bounds.min)/4);
        pivot_point[6] = floor(bounds.min + (bounds.max - bounds.min)/5);
        pivot_point[7] = floor(bounds.min + (bounds.max - bounds.min)/6);
        pivot_point[8] = floor(bounds.min + (bounds.max - bounds.min)/2);
        pivot_point[9] = floor(bounds.min + (bounds.max - bounds.min)/6);
        pivot_point[10] = floor(bounds.max - (bounds.max - bounds.min)/5);
        pivot_point[11] = floor(bounds.max - (bounds.max - bounds.min)/4);
        pivot_point[12] = floor(bounds.max - (bounds.max - bounds.min)/2);
        pivot_point[13] = floor(bounds.max - (bounds.max - bounds.min)/1.5);
        pivot_point[14] = floor(bounds.max - (bounds.max - bounds.min)/1.2);
        pivot_point[15] = floor(bounds.min);


        double A_error = pow(1,-2);
        double B_error = 1;

        double N = N_terms;

        // get the current error for the initial conditions

        //How to improve this method :
        // Check slope data, and run convergence test
        // If diverging break.

        for (int i=0;i<14;i++){
            best_error += JunctionModeler::absolute(A_best*sinh(B_best*inputData.x[int(pivot_point[i])])
                    - inputData.y[int(pivot_point[i])])/JunctionModeler::absolute(inputData.y[int(pivot_point[i])]);
        }

        // loop from 3 to 8 : range of values used to describe from where curve becomes linear to very curvy
        for (int j = 3;j<8;j++){
            // get the total error between all pivot points given current conditions
            double total_error = 0;

            // Find the current A value, and the step size for A and B
            double A_i = pow(10,-j) - A_error;
            double A_step = 2*A_error/N;
            double B_step = 2*B_error/N;

            // Loop through all A values
            for (int i=0;i<N;i++){

                // Set the current A value, such that its range matches only includes
                // regions more likely to contain the best fit conditions (formula found via excel)
                // trendline feature

                double B_i = 2*j - 6 + 0.1 ;

                // Looop through all B values in range
                for (int k=0;k<N;k++) {
                    total_error = 0;
                    // Get the total error for each B value in range
                    for (int n=0;n<14;n++){
                        total_error += JunctionModeler::absolute(A_i*sinh(B_i*inputData.x[int(pivot_point[n])])
                                - inputData.y[int(pivot_point[n])])/JunctionModeler::absolute(inputData.y[int(pivot_point[n])]);
                    }
                    // If your condition have a smaller error then previous best conditions,
                    // store the coefficients

                    if (total_error < best_error){
                        best_error = total_error;
                        A_best = A_i;
                        B_best = B_i;
                    }
                    B_i += B_step;
                }
                A_i += A_step;
            }
            A_error /= 20;
            B_error += 1;
        }

    }else{

        A_best = A_manual;
        B_best = B_manual;
        bounds.min = 0;
        bounds.max = inputData.x.size();

    }
    JunctionModeler::hyperbolic_model model;

    for (int i=0;i<inputData.x.size();i++){
        model.x.append(inputData.x[i]);
        model.y.append(A_best*sinh(B_best*inputData.x[i]));
    }


    model.A = A_best;
    model.B = B_best;
    model.minBound = bounds.min;
    model.maxBound = bounds.max;
    model.best_error = best_error;
    return model;

}


JunctionModeler::VoltagePointer amplify(double V1,double freq, double t,
                                        JunctionModeler::hyperbolic_model multiRegion[4],double step,double Rin = 10000){
    double angle = 2*3.14*t*freq;
    double next_angle = 2*3.14*(t+step)*freq;
    double current_sign = JunctionModeler::sgn(sin(angle));
    double Vin = V1*JunctionModeler::absolute(sin(angle));

    bool is_increasing = true;

    // check how the point is movign, to assign a given sinh model (ie where to get the current current ask)

    if (sin(next_angle) < sin(angle))
        is_increasing = false;


    double Itarget = Vin/Rin;

    // Check which region of the IV curve we are looking at, and calculate the best sinh fit
    // boolean's return 1 or 0, by combing two booleans ,
    // and increasing the weight of one , we can have four different cases

    int region;
    switch (is_increasing + 2*(current_sign < 0 ) ){
        case 1 :
            region = 0;
            break;
        case 0 :
            region = 2;
            break;
        case 2 :
            region = 1;
            break;
        case 3 :
            region = 3;
            break;
    }

    double Vout = asinh(Itarget/multiRegion[region].A)/multiRegion[region].B;
    if (current_sign <= 0){
        Vout *= -1;
    }

    JunctionModeler::VoltagePointer mySignal;

    mySignal.input = -1*V1*sin(2*3.14*t*freq);
    mySignal.output = Vout;

    return mySignal;
}

JunctionModeler::VectorSet get_table_from_clipboard(){
    // Get current text in clipboard
    QClipboard *clipboard = QApplication::clipboard();
    QString clipboard_text = clipboard->text();
    JunctionModeler::VectorSet set;

    // Convert data in clipboard to a vector of strings (separated by iterator \n )
    istringstream iss(clipboard_text.toStdString());
    vector<string> items;
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter(items));

    int num_items = items.size();

    // Place every element in vector into the table
    int col = 2;
    int row = -1;

    for (int i=0;i<num_items;i++){
        QString item = QString::fromStdString(items[i]);
        if ( col > 1 ){
            col = 0;
            row += 1;
            set.x.append(item.toDouble());
        }else {
            set.y.append(item.toDouble());
        }
        col +=1 ;
    }
    return set;

}
double sinh_converter(JunctionModeler::hyperbolic_model* regionOutputData, JunctionModeler::hyperbolic_model* regionInputData,
                    QVector<double> xdata, QVector<double> ydata, bool accountForHysteresis,
                    double thresholds[4],double manual_A_coeff[4],double manual_B_coeff[4], double N, bool in_manual = false, bool four_sweep_mode = true, bool no_data = false){

    double best_error = 0;

    if (!no_data) {
        JunctionModeler::VectorSet someData;
        someData.x = xdata;
        someData.y = ydata;

        JunctionModeler::SlicedVector sliced_data = slice_data(someData,0,accountForHysteresis,four_sweep_mode);
        // JunctionModeler::hyperbolic_model regionOutputData[4];
        if(sliced_data.errors) {
            regionOutputData[0].errors = true;
        }

        JunctionModeler::VectorSet vectorSets[4] = {sliced_data.upwards_positive,sliced_data.downwards_negative,
                                   sliced_data.downwards_positive,sliced_data.upwards_negative};


        if (accountForHysteresis == false){
            vectorSets[2] = vectorSets[0];
            vectorSets[3] = vectorSets[1];
            thresholds[2] = thresholds[0];
            thresholds[3] = thresholds[1];
        }

        for (int i=0;i<4;i++){
           regionInputData[i].x = vectorSets[i].x;
           regionInputData[i].y = vectorSets[i].y;

           if(regionInputData[i].x.size() < 3 ||  regionInputData[i].y.size() < 3) {
               regionOutputData[0].errors = true;
           }

           JunctionModeler::VectorSet sinh_input_data;
           sinh_input_data.x = regionInputData[i].x;
           sinh_input_data.y = regionInputData[i].y;
           if (sinh_input_data.x.size() >= 15 || sinh_input_data.y.size() >= 15){
            regionOutputData[i] = fit_sinh(sinh_input_data,thresholds[i],in_manual,manual_A_coeff[i],manual_B_coeff[i],(int)N);
            best_error += regionOutputData[i].best_error/4;
           }
        }
    }else {
        qDebug("found no data");
        for (int i=0;i<4;i++){
           regionInputData[i].x = {0,0,0};
           regionInputData[i].y = {0,0,0};

           JunctionModeler::VectorSet sinh_input_data;
           sinh_input_data.x = QVector<double>::fromStdVector(conv_to<std::vector<double>>::from(linspace(-1,1,500)));
           sinh_input_data.y = {0,0,0};
           regionOutputData[i] = fit_sinh(sinh_input_data,thresholds[i],true,manual_A_coeff[i],manual_B_coeff[i],(int)N);
        }
    }
    return best_error;

}

double generate_waveform(QVector<double>* waveform_data, double *waveform_bounds,
                                        JunctionModeler::hyperbolic_model regionOutputData[4],
                                        double freq_input, double approx_max_time, double input_resistance,
                                        double input_voltage) {

    // Find closest N-multiple of the nyquest range from the
    // input selected by the user. This minimalizes widowing effects, and makes power spectra more accurate
    JunctionModeler::VoltagePointer wavedata;

    // Using higher then 44100 would make the data very messy,
    //(higher res. != better results always when dealign with log stuff, lower FS will lower res.)
    double Fs = 44100;
    double nyquistRange = 2 * 1/freq_input;

    // Find closest N-multiple of the nyquest range from the
    // input selected by the user. This minimalizes widowing effects, and makes power spectra more accurate

    // Get closest n multiple to the nyquist range (to avoid any windowing no matter the users input)
    int closest_n_multiple = (int)round(approx_max_time/nyquistRange);
    double time_after_nyquist = nyquistRange*(double)closest_n_multiple;

    vec t = linspace(0,nyquistRange*(double)closest_n_multiple,2*Fs);

    // Sampling rate is 1 / the time difference between two samples:
    double step = t(2)-t(1);

    waveform_bounds[0] = 100000;
    waveform_bounds[1] = -100000;

    // Create table of output waveform data:
    for (int i=1;i<2*Fs;i++){
        wavedata = amplify(input_voltage,freq_input,t(i),
                           regionOutputData,step,input_resistance);
        waveform_data[0] << t(i);
        waveform_data[1] << wavedata.output;
        waveform_data[2] << wavedata.input;
        if (wavedata.output < waveform_bounds[0])
            waveform_bounds[0] = wavedata.output;
        if (wavedata.output > waveform_bounds[1])
            waveform_bounds[1] = wavedata.output;
    }
    return time_after_nyquist;
}

void getPowerSpectra(QVector<double> *power_spectra, QVector<double> waveform_data[3], double input_frequency) {

    double rate = 1/(waveform_data[0][1]-waveform_data[0][0]);

    // Take the fourier transform of the data :
    QVector<double> x_array;

    int n = waveform_data[0].size();

    // Convert the array from the data to a linear algebra defined vector :
    for (int i=0;i<n;i++)
        x_array.append(waveform_data[1][i]);

    vec X = vec(x_array.toStdVector());
    double nfft = JunctionModeler::pow2roundup(waveform_data[0].size());

    // Get the real and imaginary components of the fourier series :

    QVector<double> Y = QVector<double>::fromStdVector(
                            conv_to<std::vector<double>>::from(
                                real(fft(X, int(nfft)))));
    QVector<double> Y_img = QVector<double>::fromStdVector(
                        conv_to<std::vector<double>>::from(
                            imag(fft(X, int(nfft)))));

    // Half of the data is useless (just repeats itself)
    Y.resize(nfft/2);

    QVector<double> amplitude;
    QVector<double> freq;

    // Convert the FFT into power spectra (looks nicer) :
    // Don't plot everything , or else the peak finder will get very slow:

    double imax = (int)100*input_frequency*nfft/rate;
    for(int i=1;i< imax;i++){
         double Sff = 20*log10f(JunctionModeler::absolute(pow(Y[i],2)-pow(Y_img[i],2))/pow(nfft/2,2));
         if (Sff == Sff && Sff < 0){
            freq << i*rate/nfft;
            amplitude << Sff;
         }
    }

    power_spectra[0] = freq;
    power_spectra[1] = amplitude;
}

void JunctionModeler::on_pushButton_clicked()
{

    ui->pushButton->setText("loading ...");
    bool no_data = false;
    this->repaint();

    waveform_x_data.clear();
    waveform_y_data.clear();
    output_x_data.clear();
    output_y_data.clear();
    peaks_x_data.clear();
    peaks_y_data.clear();
    power_x_data.clear();
    power_y_data.clear();

    JunctionModeler::VectorSet input_data = get_table_from_clipboard();

    /*** 1. Copy data from clipboard into vectors for the program to use  ***/
        if (!this->useOldData) {
            // Setup table to have two columns, and a header row

            ui->tableWidget->setRowCount(1);
            ui->tableWidget->setColumnCount(2);
            // Add header to table
            ui->tableWidget->setItem(0, 0, new QTableWidgetItem("Voltage"));
            ui->tableWidget->setItem(0, 1, new QTableWidgetItem("Current"));

            ui->tableWidget->setRowCount(input_data.x.size());

            xdata.clear();
            ydata.clear();
            ui->tableWidget->clear();
            for (int i=0;i<input_data.x.size();i++){

                ui->tableWidget->setItem(i,0,new QTableWidgetItem(
                                             QString::number(input_data.x[i])));
                ui->tableWidget->setItem(i,1,new QTableWidgetItem(
                                             QString::number(input_data.y[i])));
            }
                xdata = input_data.x;
                ydata = input_data.y;
            }
        else {
            if (xdata.empty())
                no_data = true;
            this->useOldData = false;
        }


    /*** 2. Convert data into a Sinh Model(s)  ***/

        VectorSet inputData;
        inputData.x = xdata;
        inputData.y = ydata;


        double thresholds[4] = {ui->threshold_box_up_positive->value(),ui->threshold_box_down_negative->value(),
                                ui->threshold_box_down_positive->value(),ui->threshold_box_up_negative->value()};
        double manual_A_coeff[4] = {ui->A1_box->value(),ui->A2_box->value(),
                                    ui->A3_box->value(),ui->A4_box->value()};
        double manual_B_coeff[4] = {ui->B1_box->value(),ui->B2_box->value(),
                                    ui->B3_box->value(),ui->B4_box->value()};

        JunctionModeler::hyperbolic_model regionOutputData[4];
        JunctionModeler::hyperbolic_model regionInputData[4];
        double best_error = sinh_converter(regionOutputData,regionInputData,xdata,ydata,accountForHysteresis,thresholds,manual_A_coeff,
                       manual_B_coeff,ui->IV_N_box->value(),inManualMode,ui->four_sweep->isChecked(),no_data);
        ui->squared_error_label->setText(QString::number(powf(best_error,2)));
        if(regionOutputData[0].errors) {
            ui->pushButton->setText("Incorrect Data! Try Pasting Again, Or switch sweep mode. ");
            xdata.clear();
            return;
        }

        ui->IVPlot->detachItems();

        for (int i=0;i<4;i++){
           QPolygonF IVRegions,IVRegions_modeled;

           for (int j=1;j<regionOutputData[i].x.size();j++){
               if (!no_data)
                IVRegions << QPointF(regionInputData[i].x[j],regionInputData[i].y[j]);

               IVRegions_modeled << QPointF(regionOutputData[i].x[j], regionOutputData[i].y[j]);
               output_x_data.append(regionOutputData[i].x[j]);
               output_y_data.append(regionOutputData[i].y[j]);
           }
           if (!no_data){
               ui->IVPlot->replot();
               QwtPlotCurve *iv_curve = new QwtPlotCurve();
               iv_curve->setPen( Qt::blue, 4 ),
                       iv_curve->setRenderHint( QwtPlotItem::RenderAntialiased, true );
               iv_curve->setSamples(IVRegions);
               iv_curve->attach(ui->IVPlot);
           }

           ui->IVPlot->replot();
           QwtPlotCurve *iv_curve2 = new QwtPlotCurve();
           iv_curve2->setPen( Qt::red, 4 ),
                   iv_curve2->setRenderHint( QwtPlotItem::RenderAntialiased, true );
           iv_curve2->setSamples(IVRegions_modeled);
           iv_curve2->attach(ui->IVPlot);
           ui->IVPlot->replot();

        }

    /*** 3. Simulate Opamp Circuit to capture Waveforms ***/

        QPolygonF output_waveform, input_waveform;

        QVector<double> waveform_data[3];

        double waveform_bounds[2];

        // Calculate automatic R value for Amplification at Vin for testing purposes
        double auto_R1 = ui->VinBox->value()/(regionOutputData[0].A*sinh(1*regionOutputData[0].B));
        double auto_R2 = ui->VinBox->value()/(regionOutputData[1].A*sinh(1*regionOutputData[1].B));
        if (ui->Rin_auto->isChecked())
            ui->R1_box->setValue( (auto_R1 + auto_R2)/2);

        double time_after_nyquist = generate_waveform(waveform_data,waveform_bounds,regionOutputData,
                                                ui->inputFreqBox->value(), ui->maxTimeBox->value(), ui->R1_box->value(),
                                                ui->VinBox->value());
        ui->maxTimeBox->setValue(time_after_nyquist);
        // Create table of output waveform data:
        for (int i=0;i<waveform_data[0].size();i++){
            output_waveform << QPointF(waveform_data[0][i],waveform_data[1][i]);
            input_waveform << QPointF(waveform_data[0][i],waveform_data[2][i]);
        }

        // Save data as public variables used for save results in the ui
        waveform_x_data = waveform_data[0];
        waveform_y_data = waveform_data[1];


    /*** 4. Get the Power Spectra of the waveform  ***/
        QVector<double> power_spectra[2];
        QPolygonF output_fft;
        getPowerSpectra(power_spectra, waveform_data, ui->inputFreqBox->value());
        for (int i=0;i< power_spectra[0].size();i++)
            output_fft << QPointF(power_spectra[0][i],power_spectra[1][i]);

    /*** 5. Get the peaks of the power spectra  ***/
        VectorSet outputFreq = grabFrequencies(power_spectra[0], power_spectra[1],ui->fftThresholdBox->value());
        power_x_data = power_spectra[0];
        power_y_data = power_spectra[1];


    /*** 6. Add all collected data to the Graphic User Interface  ***/

        // repaint all four sets of data with their own corresponding curve
        ui->qwtPlot_2->detachItems();
        ui->qwtPlot_3->detachItems();


        QwtPlotZoomer* zoomer2 = new QwtPlotZoomer(ui->qwtPlot_2->canvas(),false);
        zoomer2->setRubberBandPen( QColor( Qt::black ) );
        zoomer2->setTrackerPen( QColor( Qt::black ) );

        QwtPlotCurve *curve3 = new QwtPlotCurve();
        curve3->setPen( Qt::red, 4 ),
                curve3->setRenderHint( QwtPlotItem::RenderAntialiased, true );
        curve3->setSamples(output_waveform);
        curve3->attach(ui->qwtPlot_2);
        ui->qwtPlot_2->replot();

        QwtPlotCurve *curve4 = new QwtPlotCurve();
        curve4->setPen( Qt::blue, 4 ),
                curve4->setRenderHint( QwtPlotItem::RenderAntialiased, true );
        curve4->setSamples(input_waveform);
        curve4->attach(ui->qwtPlot_2);
        ui->qwtPlot_2->setAxisAutoScale(ui->qwtPlot_2->yLeft);
        ui->qwtPlot_2->setAxisAutoScale(ui->qwtPlot_2->xBottom);
        ui->yLeftMin_plot2->setValue(waveform_bounds[0]);
        ui->yLeftMax_plot2->setValue(waveform_bounds[1]);
        ui->xBottomMin_plot2->setValue(waveform_data[0].first());
        ui->xBottomMax_plot2->setValue(waveform_data[0].last());
        ui->qwtPlot_2->setAxisScale(ui->qwtPlot_2->yLeft,waveform_bounds[0],waveform_bounds[1]);
        ui->qwtPlot_2->setAxisScale(ui->qwtPlot_2->xBottom,ui->xBottomMin_plot2->value(),waveform_data[0].last());
        ui->qwtPlot_2->replot();


        // Add magnification features to plot 5:
        // Allow the user to zoom in manually :
        QwtPlotZoomer* zoomer = new QwtPlotZoomer(ui->qwtPlot_3->canvas(),false);
        zoomer->setRubberBandPen( QColor( Qt::black ) );
        zoomer->setTrackerPen( QColor( Qt::black ) );


        // Setup table to have two columns, and a header row
        ui->finderTable->clear();
        ui->finderTable->setRowCount(1);
        ui->finderTable->setColumnCount(2);

        //  Add header to table
        ui->finderTable->setItem(0, 0, new QTableWidgetItem("Frequency"));
        ui->finderTable->setItem(0, 1, new QTableWidgetItem("Amplitude"));

        ui->finderTable->setRowCount(outputFreq.x.size() + 1);

        // Attatch markers using the data collected from the peak finder
        for (int i=0;i<outputFreq.x.size();i++){
            QwtPlotMarker* markers = new QwtPlotMarker();
            QwtSymbol* sym = new QwtSymbol;

            ui->finderTable->setItem(outputFreq.x.size() - i,0,new QTableWidgetItem(
                                         QString::number(outputFreq.x[i])));
            ui->finderTable->setItem(outputFreq.x.size() - i,1,new QTableWidgetItem(
                                         QString::number(outputFreq.y[i])));
            sym->setStyle(QwtSymbol::Ellipse);
            sym->setBrush(QBrush(Qt::red,Qt::SolidPattern));
            sym->setSize(10,10);
            markers->setSymbol(sym);

            markers->setValue( QPointF( outputFreq.x[i], outputFreq.y[i] ) );
            markers->attach(ui->qwtPlot_3 );
        }
        peaks_x_data  = outputFreq.x;
        peaks_y_data = outputFreq.y;

        QwtPlotCurve *curve5 = new QwtPlotCurve();
        curve5->setPen( Qt::blue, 4 ),
                curve4->setRenderHint( QwtPlotItem::RenderAntialiased, true );
        curve5->setSamples(output_fft);
        curve5->attach(ui->qwtPlot_3);
        ui->qwtPlot_3->setAxisAutoScale(ui->qwtPlot_3->yLeft);
        ui->qwtPlot_3->setAxisAutoScale(ui->qwtPlot_3->xBottom);
        ui->qwtPlot_3->replot();
        ui->yLeftMin_plot3->setValue(outputFreq.y[0]);
        ui->xBottomMax_plot3->setValue(outputFreq.x[0]);
        ui->qwtPlot_3->setAxisScale(ui->qwtPlot_3->yLeft,outputFreq.y[0],ui->yLeftMax_plot3->value());
        ui->qwtPlot_3->setAxisScale(ui->qwtPlot_3->xBottom,ui->xBottomMin_plot3->value(),outputFreq.x_max);

        ui->pushButton->setText("Paste from clipboard");

        // display coefficients for sinhmodels
        ui->A1_box->setValue(regionOutputData[0].A);
        ui->B1_box->setValue(regionOutputData[0].B);
        ui->A2_box->setValue(regionOutputData[1].A);
        ui->B2_box->setValue(regionOutputData[1].B);
        ui->A3_box->setValue(regionOutputData[2].A);
        ui->B3_box->setValue(regionOutputData[2].B);
        ui->A4_box->setValue(regionOutputData[3].A);
        ui->B4_box->setValue(regionOutputData[3].B);
}

void JunctionModeler::on_checkBox_clicked(bool checked)
{
    ui->hysteresis_group->setEnabled(checked);
    this->accountForHysteresis = checked;
}

void JunctionModeler::on_checkBox_2_clicked(bool checked)
{
    ui->manualFunction_group->setEnabled(checked);
    this->inManualMode = checked;
}

void JunctionModeler::on_rescale_power_clicked()
{
    ui->qwtPlot_3->setAxisScale(ui->qwtPlot_3->yLeft,ui->yLeftMin_plot3->value(),ui->yLeftMax_plot3->value());
    ui->qwtPlot_3->setAxisScale(ui->qwtPlot_3->xBottom,ui->xBottomMin_plot3->value(),ui->xBottomMax_plot3->value());
}

void JunctionModeler::on_rescale_power2_clicked()
{
    ui->rescale_power->click();
}

void JunctionModeler::on_exportIV_clicked()
{
    exportToDesktop("IVCurve",output_x_data,output_y_data,ui->dataPrefixInput->text(),ui->folderCheckbox->isChecked());
}

void JunctionModeler::on_exportPeaks_clicked()
{
    exportToDesktop("Peaks",peaks_x_data,peaks_y_data,ui->dataPrefixInput->text(),ui->folderCheckbox->isChecked());
}

void JunctionModeler::on_exportWaveform_clicked()
{
     exportToDesktop("Waveform",waveform_x_data,waveform_y_data,ui->dataPrefixInput->text(),ui->folderCheckbox->isChecked());
}

void JunctionModeler::on_exportPower_clicked()
{
    exportToDesktop("PowerSpectra",power_x_data,power_y_data,ui->dataPrefixInput->text(),ui->folderCheckbox->isChecked());

}

void JunctionModeler::on_rescale_waveform_clicked()
{
    ui->qwtPlot_2->setAxisScale(ui->qwtPlot_2->yLeft,ui->yLeftMin_plot2->value(),ui->yLeftMax_plot2->value());
    ui->qwtPlot_2->setAxisScale(ui->qwtPlot_2->xBottom,ui->xBottomMin_plot2->value(),ui->xBottomMax_plot2->value());
    ui->qwtPlot_2->replot();
}

void JunctionModeler::on_rescale_waveform2_clicked()
{
    ui->rescale_waveform->click();
}

void JunctionModeler::on_exportSettings_clicked()
{
    QPixmap screenshot = QPixmap::grabWidget(window());
    QString fname;
    QString prefix = ui->dataPrefixInput->text();
    bool inFolder = ui->folderCheckbox->isChecked();
    if (inFolder){
        if (!QDir(QStandardPaths::writableLocation(QStandardPaths::DesktopLocation) + "\\" +prefix).exists())
            QDir().mkdir(QStandardPaths::writableLocation(QStandardPaths::DesktopLocation) + "\\" +prefix);
        fname = QStandardPaths::writableLocation(QStandardPaths::DesktopLocation) + "\\" +prefix + "\\"+ prefix + "_screenshot.png";
    }else {
         fname = QStandardPaths::writableLocation(QStandardPaths::DesktopLocation) + "\\"+ prefix + "_screenshot.png";
    }

    screenshot.save(fname);


}

void JunctionModeler::tutorial() {
    //Get current text in clipboard
    char temp;

    cout << "Junctions become linear at high enough voltages." << endl;
    cout << "Enter a threshold of confidence required to consider data to be linear " << endl;
    cout << "threshold increases logarithmically with range 0.0000001 to 50." << endl;
    cout << "Input a value within this range (recommend 0.001) " << endl;
    double threshold;
    cin >> threshold;
    cout << "The threshold you entered is : "<< threshold << endl;
    cout << "If you find that your data is inaccurate at high or low voltages " << endl;
    cout << ", you can adjust this in a future run. Also you can enter multiple thresholds to account for asymmetry." << endl;
    cout << "::: Get IV Curve :::" << endl;
    cout << "Please copy IV curve from excel to clipboard, once completed" << endl;
    cout << "press y to confirm (you do not need to paste! )" << endl;

    cin >> temp;

    /*** 1. Grab data from clipboard  ***/
    JunctionModeler::VectorSet set = get_table_from_clipboard();
    exportToDesktop("IVCurve",set.x,set.y,QString::fromStdString(argv[3]),true,
            QString::fromStdString(argv[2]));

    double thresholds[4] = {threshold,threshold,
                            threshold,threshold};
    double manual_A_coeff[4] = {0,0,
                               0,0};
    double manual_B_coeff[4] = {0,0,
                                0,0};

    JunctionModeler::hyperbolic_model regionOutputData[4];
    JunctionModeler::hyperbolic_model regionInputData[4];
    double best_error = sinh_converter(regionOutputData,regionInputData,set.x,set.y,false,thresholds,manual_A_coeff,
                   manual_B_coeff,280,false);
    JunctionModeler::VectorSet iv_out;
    for (int i=0;i<4;i++){
        for (int j=1;j<regionOutputData[i].x.size();j++){
            iv_out.x.append(regionOutputData[i].x[j]);
            iv_out.y.append(regionOutputData[i].y[j]);
        }
    }

    /*** 2. Calculate the fit for the IV Curve  ***/
    exportToDesktop("IVCurve_fitted",iv_out.x,iv_out.y,QString::fromStdString(argv[3]),true,
            QString::fromStdString(argv[2]));
    cout << "Please confirm that the data and the fitted data were exported successfully" << endl;
    cout << "You may want to check the fit accuracy in excel before continuing." << endl;
    cout << "To do this plot both sets of data ontop of eachother using Scatter with only Markers" << endl;
    cout << "If you are confident with the fit provided , press y to continue" << endl;
    cin >> temp;
    cout << "::: Creating a Waveform :::" << endl;
    cout << "The program will now put the simulated device through an opamp circuit." << endl;
    cout << "As you are in guided/easy mode we will assume that you are sending an input signal of :" << endl;
    cout << "0.3V, and amplifying the signal to 1V. If you passed in Vin and Rin during execution" << endl;
    cout << "this will not be the case. Furthermore it is assumed you are using an input frequency of " << endl;
    cout << "150Hz, and are sampling over a range of 0.226667 seconds" << endl;

    QVector<double> waveform_data[3];
    double waveform_bounds[2];

    // Calculate automatic R value for Amplification at Vin for testing purposes
    double auto_R1 = 0.3/(regionOutputData[0].A*sinh(1*regionOutputData[0].B));
    double auto_R2 = 0.3/(regionOutputData[1].A*sinh(1*regionOutputData[1].B));
    double auto_R = (auto_R1 + auto_R2)/2;

    double time_after_nyquist = generate_waveform(waveform_data,waveform_bounds,regionOutputData,
                                            150, 0.226667, auto_R, 0.3);

    /*** 3. Get the Waveform of the IV data  ***/
    exportToDesktop("Waveform",waveform_data[0],waveform_data[1],QString::fromStdString(argv[3]),true,
            QString::fromStdString(argv[2]));
    cout << "The waveform was exported to deskop, to acheive a 1V range a resistance of : " << endl;
    cout << auto_R << " was inputed into the opamp circuit." << endl;
    cout << "The max and min of the waveform are " << waveform_bounds[0] << " and " << waveform_bounds[1] << endl;
    cout << "Note that these will not be 1V exactly if your junction is not symmetric." << endl;
    cout << "The time range was adjusted to : " << time_after_nyquist << " to avoid windowing effects. " << endl;
    cout << "If you are happy with the waveform exported, press y to continue" << endl;
    cin >> temp;
    cout << ":::Capturing Power Spectra:::" << endl;
    cout << "Caluating the variance of the data distributed over the frequency components." << endl;

    /*** 4. Get the Power Spectra of the waveform  ***/
    QVector<double> power_spectra[2];
    getPowerSpectra(power_spectra, waveform_data, 150);
    exportToDesktop("Power Spectra",power_spectra[0],power_spectra[1],QString::fromStdString(argv[3]),true,
            QString::fromStdString(argv[2]));

    /*** 5. Get the peaks of the power spectra  ***/
    VectorSet outputFreq = grabFrequencies(power_spectra[0], power_spectra[1],50);

    cout << "The peaks of the power spectra will be identified using the data collected." << endl;
    cout << "A threshold value is used to determine the requirement for a maximum to be" << endl;
    cout << "considered a peak, this value ranges from 0 to 100. " << endl;
    cout << "In guided mode this is set automatically to : 50" << endl;

    exportToDesktop("Spectra Peaks",outputFreq.x,outputFreq.y,QString::fromStdString(argv[3]),true,
            QString::fromStdString(argv[2]));

    cout << "All data should be generated!" << endl;
    cout << "Thank you for following the guided tour of the JunctionModeler Library" << endl;
    cout << "All functions were developed by Lucas Zeer and are property of NINT. " << endl;
    cout << "When ready, feel free to close the command window." << endl ;

}

void JunctionModeler::on_exportAll_clicked()
{
    ui->exportSettings->click();
    ui->exportIV->click();
    ui->exportWaveform->click();
    ui->exportPower->click();
    ui->exportPeaks->click();
}

void JunctionModeler::on_useOldData_clicked()
{
    this->useOldData = true;
    ui->pushButton->click();
}

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <iterator>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "cluster.h"

#define PY_SSIZE_T_CLEAN


using namespace std;


//A CLUSTER OBJECT ("cluster"):
Cluster cluster;



/**
   *  Get the energy of the cluster with the greatest energy. - 
   * 
   *  [Compares energies of all clusters from layer 1 and returns the greatest sum]
   * 
   *  @return  (double) the sum of the energies of all pixels in that cluster.
  */
double get_maximum_total_energy()
{
    double eMax = 0;
    double j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            eMax = (cluster.get_total_energy(i));
        }
        else if ((j = cluster.get_total_energy(i)) > eMax)
        {
            eMax = j;
        }
    }

    return eMax;
}

/**
   *  Get the energy of the cluster with the smallest energy. - 
   * 
   *  [Compares energies of all clusters from layer 1 and returns the smallest sum]
   * 
   *  @return  (double) the sum of the energies of all pixels in that cluster.
  */
double get_minimum_total_energy()
{
    double eMax = 0;
    double j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            eMax = (cluster.get_total_energy(i));
        }
        else if ((j = cluster.get_total_energy(i)) < eMax)
        {
            eMax = j;
        }
    }

    return eMax;
}

/**
   *  Get the delta time of the cluster with the smallest delta time. - 
   * 
   *  [Compares delta times of all clusters from layer 1 and returns the smallest one]
   * 
   *  @return  (double) the smallest delta time.
  */
double get_minimum_time()
{
    double tMin = 0;
    double j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            tMin = cluster.get_delta_time(i);
        }
        else if ((j = cluster.get_delta_time(i)) < tMin)
        {
            tMin = j;
        }
    }

    return tMin;
}

/**
   *  Get the size of the cluster with the greatest size. - 
   * 
   *  [Compares sizes of all clusters from layer 1 and returns the greatest one]
   * 
   *  @return  (int) the size of the cluster with the greatest size.
  */
int get_maximum_size()
{
    int sMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            sMax = (cluster.get_size(i));
        }
        else if ((j = cluster.get_size(i)) > sMax)
        {
            sMax = j;
        }
    }

    return sMax;
}

/**
   *  Get the size of the cluster with the smallest size. - 
   * 
   *  [Compares sizes of all clusters from layer 1 and returns the smallest one]
   * 
   *  @return  (int) the size of the cluster with the smallest size.
  */
int get_minimum_size()
{
    int sMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            sMax = (cluster.get_size(i));
        }
        else if ((j = cluster.get_size(i)) < sMax)
        {
            sMax = j;
        }
    }

    return sMax;
}

/**
   *  Get the delta time of the cluster with the greatest delta time. - 
   * 
   *  [Compares delta times of all clusters from layer 1 and returns the greatest one]
   * 
   *  @return  (double) the greatest delta time.
  */
double get_maximum_time()
{
    double tMin = 0;
    double j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            tMin = cluster.get_delta_time(i);
        }
        else if ((j = cluster.get_delta_time(i)) > tMin)
        {
            tMin = j;
        }
    }

    return tMin;
}

/**
   *  Get the x-coordinate size difference of the cluster with the smallest x-coordinate size difference. - 
   * 
   *  [Compares x-coordinate size differences of all clusters from layer 1 and returns the smallest one]
   * 
   *  @return  (double) the smallest x-coordinate size difference.
  */
int get_minimum_x_dimension()
{
    int dMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            dMax = (cluster.get_x_dimension(i));
        }
        else if ((j = cluster.get_x_dimension(i)) < dMax)
        {
            dMax = j;
        }
    }

    return dMax;
}

/**
   *  Get the x-coordinate size difference of the cluster with the greatest x-coordinate size difference. - 
   * 
   *  [Compares x-coordinate size differences of all clusters from layer 1 and returns the greatest one]
   * 
   *  @return  (double) the greatest x-coordinate size difference.
  */
int get_maximum_x_dimension()
{
    int dMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            dMax = (cluster.get_x_dimension(i));
        }
        else if ((j = cluster.get_x_dimension(i)) > dMax)
        {
            dMax = j;
        }
    }

    return dMax;
}

/**
   *  Get the y-coordinate size difference of the cluster with the smallest y-coordinate size difference. - 
   * 
   *  [Compares y-coordinate size differences of all clusters from layer 1 and returns the smallest one]
   * 
   *  @return  (double) the smallest y-coordinate size difference.
  */
int get_minimum_y_dimension()
{
    int dMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            dMax = (cluster.get_y_dimension(i));
        }
        else if ((j = cluster.get_y_dimension(i)) < dMax)
        {
            dMax = j;
        }
    }

    return dMax;
}

/**
   *  Get the y-coordinate size difference of the cluster with the greatest y-coordinate size difference. - 
   * 
   *  [Compares y-coordinate size differences of all clusters from layer 1 and returns the greatest one]
   * 
   *  @return  (double) the greatest y-coordinate size difference.
  */
int get_maximum_y_dimension()
{
    int dMax = 0;
    int j = 0;

    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        if (i == 0)
        {
            dMax = (cluster.get_y_dimension(i));
        }
        else if ((j = cluster.get_y_dimension(i)) > dMax)
        {
            dMax = j;
        }
    }

    return dMax;
}


//IN ORDER TO MAKE THIS CODE WORK, IT IS NECESSARY TO FILL THE ABSOLUTE PATHS IN MAIN FOR "LOAD_DATA, OUTPUT AND EXPORT"!


//MAIN  
int main(int arg, char *argv[])
{
    // index of the cluster in the vector of clusters [vector<vector<Pixel*>> clusters]
    int n = 0;
    // the pixel atributes to be read from the source file [x-coordinate, y-coordinate, time, energy]
    std::vector<bool> params = {true, true, true, true};


    //LOAD CLUSTERED DATA INTO LAYER 1 CLUSTERS
    cluster.l1_clusters = cluster.load_data("HERE FILL THE ABSOLUTE PATH", 5000, params, false);


    //ASK FOR THE INDEX OF THE CLUSTER IN LAYER 1 TO BE EXAMINED
    cout << "Write the number of the cluster you want to examine: \t";
    cin >> n;


    //ANALYZE A GIVEN CLUSTER AND PRINT SOME OF ITS ATRIBUTES
    printf("\n total energy: \t %f", cluster.get_total_energy(n));
    printf("\n minimum time: \t %lf", cluster.get_minimum_time(n));
    printf("\n x-coordinate difference in size: \t %i", cluster.get_x_dimension(n));
    printf("\n y-coordinate difference in size: \t %i", cluster.get_y_dimension(n));
    printf("\n x-coordinate of pixel with minimum time: \t %i", cluster.get_pixel_with_minimum_time(n).x);
    printf("\n y-coordinate of pixel with minimum time: \t %i", cluster.get_pixel_with_minimum_time(n).y);
    printf("\n x-coordinate of pixel with maximum time: \t %i", cluster.get_pixel_with_maximum_time(n).x);
    printf("\n y-coordinate of pixel with maximum time: \t %i", cluster.get_pixel_with_maximum_time(n).y);
    printf("\n delta time: \t %lf", cluster.get_delta_time(n));
    cout << endl;
    cout << endl;


    //ANALYZE THE WHOLE FILES AND PRINT SOME OF THE ATRIBUTES IN A PSEUDO-TABLE
    printf("\t Attribute \t \t Min \t \t \t Max \t");
    printf("\n Size [in pixels] \t %i \t \t \t %i \t \t", get_minimum_size(), get_maximum_size());
    printf("\n Total energy \t \t %lf \t \t %lf \t \t", get_minimum_total_energy(), get_maximum_total_energy());
    printf("\n Delta time \t \t %lf \t \t %lf \t \t", get_minimum_time(), get_maximum_time());
    printf("\n X dimension \t \t %i \t \t \t %i \t \t", get_minimum_x_dimension(), get_maximum_x_dimension());
    printf("\n Y dimension \t \t %i \t \t \t %i \t \t", get_minimum_y_dimension(), get_maximum_y_dimension());

    cout << endl;
    cout << endl;


    //OUTPUT GIVEN DATA INTO AN OUTPUT FILE
    cluster.output("HERE FILL THE ABSOLUTE PATH", "get_dimensions");


    //EXPORT DATA FOR MS EXCEL ANALYSIS
    cluster.exp("HERE FILL THE ABSOLUTE PATH");


    //PRINT LINEARITY OF FIRST 2O CLUSTER FORM LAYER 1
    for (int i = 0; i < 20; i++)
    {
        printf("Cluster %i: \t LINEARITY: %f \t ROUNDNESS: %f \n", i, cluster.get_linearity(i), cluster.get_roundness(i));
    }
    cout << endl;
    cout << endl;


    //PRINT THE INDICIES ALL THE LAYER 1 - LAYER 2 PAIRS OF CLUSTERS (very very slow algorithm!)
    int m = 0;  //a variable for counting the number of items that were printed
    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        cluster.get_pairs(i);
        if (!cluster.pairs[i].empty())
        {
            for (int j = 0; j < cluster.pairs[i].size(); j++)
            {
                printf("Pair #%i: \t l1_clusters[%i], \t l2_clusters[%i] \n", m, i, cluster.pairs[i][j]);
                m++;
            }
        }
    }
    cout << endl;
    cout << endl;


    //PRINT THE INDICIES OF THE PIXELS WITH THE SMALLEST TIMES FROM ALL THE CLUSTER FROM LAYER 1
    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        cluster.get_minimum_time_index(i, "l1_clusters");
    }
    cout << endl;
    cout << endl;


    //PRINT THE DISTANCE BETWEEN THE TWO DETECTOR LAYERS USING DIFFERENT TECHNIQUES AND DIFFERENT CLUSTERS
    int o = 0;  //a variable used find two consecutive pairs for the analysis of the distance between the two detector layers
    vector<int> indicies;   //a vector of all indicies of l1_clusters that are in a pair with another cluster form l2_clusters
    float distance; //the distance between the two detector layers
    float height;   //the height between the two detector layers
    m = 0;
    for (int i = 0; i < cluster.l1_clusters.size(); i++)
    {
        cluster.get_pairs(i);
        if (!cluster.pairs[i].empty())
        {
            for (int j = 0; j < cluster.pairs[i].size(); j++)
            {
                printf("pair #%i (%i) - is pair: \t %s \n", i, m, "true");
                indicies.push_back(i);
                if (o >= 1)
                {
                    distance = cluster.get_distance(i);
                    printf("cluster #%i: \tDISTANCE is %f\n", i, distance);
                    height = cluster.get_height(indicies[o - 1], indicies[o]);
                    printf("cluster #%i: \tHEIGHT is %f\n", i, height);
                }
                m++;
                o++;
            }
        }
    }
    cout << endl;
    cout << endl;

    return 1;
}
#include <udf.h>
#include <sys/stat.h>

/* (C)  2022
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

/*************************************************************************************/
DEFINE_ON_DEMAND(CFD_convert_csv2ip)
/*************************************************************************************/
{
#if RP_HOST

#if 1   /* local macros */

#define     _csv_file_number_start      1000
#define     _csv_file_number_delta      1000
#define     _csv_file_number_of_files   99
#define     _csv_file_number_of_lines   80000

#endif

#if 1   /* local vars */

    short   csv_file_exist;

    int     csv_file_number, ip_file_number, i_data;

    double  x[3], alpha_gas, u_gas[3], u_solid[3];

    double  ip_x[_csv_file_number_of_lines][3], ip_alpha_gas[_csv_file_number_of_lines], ip_u_gas[_csv_file_number_of_lines][3], ip_u_solid[_csv_file_number_of_lines][3];

    char    csv_filename[40], ip_filename[40];

    FILE    *f_csv = NULL, *f_ip = NULL;

#endif

    /* create output directory */
    mkdir("./data/ip",0777);

    csv_file_number = _csv_file_number_start;

    ip_file_number = 0;

    csv_file_exist = 1;

    while(csv_file_exist){

        sprintf(csv_filename, "./data/csv/%d.txt", csv_file_number);

        f_csv = fopen(csv_filename,"r");

        if(f_csv){

            /* read csv data, close f_csv, increase csv_file_number */
            {
                for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                    fscanf(f_csv, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",

                        &x[0], &x[1], &x[2], &alpha_gas, &u_gas[0], &u_gas[1],  &u_gas[2], &u_solid[0], &u_solid[1],  &u_solid[2]);

                        ip_x[i_data][0] =   x[0];
                        ip_x[i_data][1] =   x[1];
                        ip_x[i_data][2] =   x[2];

                        ip_alpha_gas[i_data] = alpha_gas;

                        ip_u_gas[i_data][0] =   u_gas[0];
                        ip_u_gas[i_data][1] =   u_gas[1];
                        ip_u_gas[i_data][2] =   u_gas[2];

                        ip_u_solid[i_data][0] =     u_solid[0];
                        ip_u_solid[i_data][1] =     u_solid[1];
                        ip_u_solid[i_data][2] =     u_solid[2];
                }

                fclose(f_csv);

                csv_file_number += _csv_file_number_delta;
            }

            /* open ip file and write ip header & data */
            {
                sprintf(ip_filename,"./data/ip/%04d.ip", ip_file_number);

                f_ip = fopen(ip_filename,"w");

                if(f_ip){

                    /* write ip header */
                    {
                        fprintf(f_ip, "3\n");
                        fprintf(f_ip, "3\n");
                        fprintf(f_ip, "%d\n", _csv_file_number_of_lines);
                        fprintf(f_ip, "8\n");
                        fprintf(f_ip, "mp-1\n");
                        fprintf(f_ip, "mp-2\n");
                        fprintf(f_ip, "x-velocity-1\n");
                        fprintf(f_ip, "y-velocity-1\n");
                        fprintf(f_ip, "z-velocity-1\n");
                        fprintf(f_ip, "x-velocity-2\n");
                        fprintf(f_ip, "y-velocity-2\n");
                        fprintf(f_ip, "z-velocity-2\n");
                    }

                    /* write ip_data */
                    {
                        /* x-coord */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_x[i_data][0]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* y-coord */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_x[i_data][1]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* z-coord */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_x[i_data][2]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* mp-1 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_alpha_gas[i_data]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* mp-2 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", (1.0 - ip_alpha_gas[i_data]));
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* x_velocity-1 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_gas[i_data][0]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* y_velocity-1 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_gas[i_data][1]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* z_velocity-1 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_gas[i_data][2]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* x_velocity-2 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_solid[i_data][0]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* y_velocity-2 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_solid[i_data][1]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                        /* z_velocity-2 */
                        {
                            fprintf(f_ip, "(");

                            for(i_data = 0; i_data < _csv_file_number_of_lines; i_data++){

                                fprintf(f_ip, "%f\n", ip_u_solid[i_data][2]);
                            }

                            fprintf(f_ip, ")\n");
                        }

                    }

                    fclose(f_ip);

                    ip_file_number++;
                }
            }
        }
        else{
            csv_file_exist = 0;
        }

        if(csv_file_number == (_csv_file_number_start + _csv_file_number_of_files * _csv_file_number_delta)){

            csv_file_exist = 0;
        }
    }

    Message("\n\n...CFD_convert_csv2ip: converted %d files\n", ip_file_number);

#endif
}


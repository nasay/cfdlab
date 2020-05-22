#ifndef _CHECKS_H_
#define _CHECKS_H_

/*  Checks if each cell in colideField
    *** If the density is betwen 0.9 and 1.1
    *** If the norm of velocity is not greater than the speed of sound
    *lenght is a vector with the size of the field and t is the current timestep
*/
void check_in_rank(Fields &fields, int *length, int t);

/*  Checks if for each cell in flagField with flag flag1 is adjacent
    *(in one of the 19 lattices) with a cel of flag flag2
    *lenght is a vector with the size of the field and t is the current timestep
*/
void check_flags(Fields &fields, int *length, int flag1, int flag2, int t);

/*  Prints out the total mass in the cavity using the data in massField
    *lenght is a vector with the size of the field and t is the current timestep
*/
void check_mass(Fields &fields, int *length, int t);

/*  Carrys al the checks (check_in_rank, check_flags, check_mass) */
void run_checks(Fields &fields, int *length, int t);

#endif

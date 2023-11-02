# t-issm: Glacier fun with J-man and Felicity 
Email felicity.holmes@geo.su.se and jamie.barnett@geo.su.se if you have questions

# Task 1: Choose your pokemon (glacier)

- For this exercise, each group will choose one of the following glaciers to model:
  - 79 North glacier (Nioghalvfjerdsbrae)
  - Helheim glacier
  - Kangerlussuaq glacier
  - Petermann glacier
  - Jakobshavn Isbrae
  - Tracy and Heilprin glaciers
 
- Each group should choose a different glacier - and its first come first served, so let Jamie + Felicity know as soon as you have chosen your fave
- Ryder glacier is only used for the in-class example, and is shown in all the example figures in these instructions

# Task 2: Run model to spun-up state 

- Open up matlab, and navigate to the t-issm directory and open the 'runme.m' file - this is where you will work from for this exercise
- Make sure all the ISSM source code AND the t-issm directory is on your Matlab PATH

## Step 1: Mesh it

- Meshing the glacier is 'Step 1' in the runme script - to run this step, you should make sure the first line of the script looks like this:

`steps = [1];`

- Once this is set, you can run the script by either pressing the large 'Run' button at the top of the screen, or by executing 'runme' in the command window

- In this step, we create a mesh from the outline provided in the accompanying shapefiles
- The minimum and maximum mesh resolutions are pre-set for each glacier, so you do not need to change these
- We use present day observed surface velocities to determine where the maximum mesh resolution is applied - with a higher mesh resolution being used in areas with faster flow

- Upon successful completion of this step, you should get output that looks similar to the example below:

<img width="1582" alt="Screenshot 2023-11-02 at 14 42 34" src="https://github.com/jamiewbarnett/t-issm/assets/141425558/d4eff310-1459-4eef-be7b-5b4d19758579">

## Step 2: Parameterise it

- In the parameterisation step, we use various observational datasets to set the **initial conditions** and **boundary conditions** for the glacier
- All of these datasets are read in and applied to the glacier using the `Greenland.par` file
- Upon successful completion of this step, you should get output that looks similar to the example below:

  <img width="1582" alt="Screenshot 2023-11-02 at 14 48 41" src="https://github.com/jamiewbarnett/t-issm/assets/141425558/2d945cbe-7184-42a1-9198-09dc2b30b408">


## Step 3: Invert it/ run the stress balance

## Step 4: Spin it up








glacier go bye-bye 

# t-issm: Glacier fun with J-man and Felis 
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

- Set `steps = [2];`
- In the parameterisation step, we use various observational datasets to set the **initial conditions** and **boundary conditions** for the glacier
- All of these datasets are read in and applied to the glacier using the `Greenland.par` file
- Upon successful completion of this step, you should get output that looks similar to the example below:

  <img width="1582" alt="Screenshot 2023-11-02 at 14 48 41" src="https://github.com/jamiewbarnett/t-issm/assets/141425558/2d945cbe-7184-42a1-9198-09dc2b30b408">

## Step 3: Invert it/ run the stress balance

- Set `steps = [3];`
- The next step is to use observed surface velocities to invert for basal friction, whilst also solving for a Stressbalance solution of the glacier
- Upon successful completion of this step, you should get output that looks similar to the example below:

 <img width="1582" alt="Screenshot 2023-11-06 at 15 02 43" src="https://github.com/jamiewbarnett/t-issm/assets/141425558/f4af377f-8a03-4e1a-8f0f-d33f4cc72675">



## Step 4: Spin it up


# Task 3: Design and run some transient simulations

- For the transient simulations, you need to decide as a group what you want to investigate
- Examples of questions you could answer are things like:
  - What are the differences in glacier behaviour up to 2050 under a low vs high emissions future scenario?
  - Is the glacier more sensitive to changes in ocean thermal forcing or to changes in SMB?
  - Under which climate scenarios does the glacier lose its floating ice tongue?
 
- You need to think about which forcings you can change, and look at the datasets available
- Some things you can change are SMB (e.g. different SSPs), basal melt under ice tongues, frontal melt along grounded termini, calving stress threshold
- It is also a good idea to do a bit of reading about your chosen glacier to understand which type of question might be relevant and/or interesting

- Once you have decided on the question(s) to investigate, you need to think about what forcings/parameters you need to change in the `runme` script, and how many simulations you want to run

- Make sure to keep good notation about which simulations you run 

# Task 4: Analyse your results

# Task 5: Make a presentation

- Each group should make a c. 20-25 minute presentation about their modelling project
- Include the following in your presentations:
    - Background on your glacier
    - Scientific questions
    - Model set-up xperimental design; how did your set up your model to address your chosen question(s)?
      - Also include information on forcings here (e.g. SMB, melt rates)
    - Results
    - Conclusions and what you would do next if you had more time

- We will then have 5-10 minutes of questions after each presentation
 

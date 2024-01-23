# t-ISSM: Glacier fun with J-man and Felis 
Email jman.barnett@geo.su.se and felicity.holmes@geo.su.se if you have questions

# Task 1: Choose your pokemon (glacier)

- For this exercise, each group will choose one of the following glaciers to model:
  - 79 North glacier (Nioghalvfjerdsbrae)
  - Helheim glacier
  - Kangerlussuaq glacier
  - Petermann glacier
  - Jakobshavn Isbrae
  - Tracy and Heilprin glaciers
 
- Each group should choose a different glacier - this will be decided in class 
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

<img width="1582" alt="Screenshot 2023-11-02 at 14 42 34" src="./Screenshots/step1.png">

## Step 2: Parameterise it

- Set `steps = [2];`
- In the parameterisation step, we use various observational datasets to set the **initial conditions** and **boundary conditions** for the glacier
- All of these datasets are read in and applied to the glacier using the `Greenland.par` file
- Upon successful completion of this step, you should get output that looks similar to the example below:

  <img width="1582" alt="Screenshot 2023-11-02 at 14 48 41" src="./Screenshots/step2.png">

## Step 3: Invert it/ run the stress balance

- Set `steps = [3];`
- The next step is to use observed surface velocities to invert for basal friction, whilst also solving for a Stressbalance solution of the glacier
- Upon successful completion of this step, you should get output that looks similar to the example below:

 <img width="1582" src="./Screenshots/step3.png">

## Step 4: Spin it up

- In the spin-up stage, we try to set parameters such as those relating to calving and frontal melt at values which lead to a stable glacier
- This stable position is then our starting point for any future transient simulations
- Spinning up a glacier is useful as it means that the way your glacier behaves in the model once you implement your forcings can be likely attributed to said forcings rather than to model drift
- The spin-up parameters for each glacier have been set for you, and these sounds give acceptable spun-up glaciers
- If this was to be a 'real' study, you would probably spend a lot more time on the spin-up and could no doubt find a set of parameters that yield an even better spun-up state

- At the end of this step, you should see something like the following:

# Task 3: Design and run some transient simulations

- For the transient simulations, you need to decide as a group what you want to investigate
- Examples of questions you could answer are things like:
  - What are the differences in glacier behaviour up to 2050 under a low vs high emissions future scenario?
  - Is the glacier more sensitive to changes in ocean thermal forcing or to changes in SMB?
  - Under which climate scenarios does the glacier lose its floating ice tongue?
 
- You need to think about which forcings you can change, and look at the datasets available
- Some things you can change are:
    - SMB (e.g. different SSPs)
    - Basal melt under ice tongues
    - Frontal melt along grounded termini
    - Calving stress thresholds for floating and grounded ice
  
- It is a good idea to do a bit of reading about your chosen glacier to understand which type of question might be relevant and/or interesting
- You also need to consider how many simulations you need to run to answer your research question

- Once you have decided on the question(s) to investigate, you need to change the relevant parameters in the header section of the `runme` script

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
 

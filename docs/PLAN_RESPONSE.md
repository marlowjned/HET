# Response to initial plan

Here I will only annotate changes I believe need to be made, or elements that need to be discussed.

# 4.3 Materials
No need for MATERIALS.md, just make the material library something easily readable (i.e. YAML, the goal should be to allow materials to define any number of properties, and trust the users to label the properties correctly, perhaps we could consider name aliasing...). This seems to me to be a poorly fleshed out section. We need to decide how this data is stored.

# 4.5 Inheritance
Where is the parent data stored? Do we have sub-configuraiton files per part? Does this design config link them all together then? Where do we decide what parts should have their own configs, and what parts should be somewhat defined within the design-level config? This section seems weak to me, and easy to get disorganized/make configs that aren't broadly useful.

# 4.6 CAD
The team is moving away from using Onshape and toward Fusion. We'll have to rely more on scripted utilities for geometry sweeps. This likely deserves some of it's own thought, and can really be a whole other sub-project, figuring out how to best sweep/generate geometries (geoemtries generated based on fixed config vars + changing runsheet params?).

# 5.2 Study vars 
This section is somewhat confusing to me. Where are categorical definitions made? Is this at the component level? Where does 'scale' become useful? What is defining the number of swept values per variable? 

# 5.3 Runsheet
Explain parameter vs study-variable space. When would it be useful for constraints to change per run? If they do change, why not define a parameter that affects how such constraints come into play (shifts the constraining conditions)? 

# 6 Versioning
Do not make versions integers, make them strings (i.e. "1.0.0", "6.2.3", etc). This version number could be printed on startup of "hetsim". I like the name "hetsim" for now, but I'll ponder a better name. Remind me before we execute.

# 7
I need more clarification on this section overall.

# 10.1 Study spec
Again, I'm suspicious of this organization, but I like the idea overall. 

# 10.3 Driver loop
We should also consider parallelization logic.

# 11 Execution
I'll trust you here, this isn't my area of expertise.

# 16-END
Much of this changes slightly with our new decisions.

# 19 Your questions
1) I like z+ as the exit direction.
2) Fixed turns, though this isn't all that important right now.
3) If the user can define a component fully, let them. Either should be possible. 
4) Well, this really comes down to how we treat the probed data. These probe selections should be easy to add/modify, so ideally make them associated with their own scripts. I.e. area selections could in theory let you also process data at the boundaries.
5) WSL is good
6) Not important right now. We can clean up the repo later.
7) Let me ponder this still, remind me to answer.


OVERALL, I like the central ideas, but there are just certain points that could be much more modular/robust. We want it to be possible to copy project data over to new projects easily, and not have to rewrite too much on deisgning a new trade, defining a new probe, etc. Let's discuss the plan further. 
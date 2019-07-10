import visual
import space_python2 as s

#define field
vec=s.vector(3,3,3)
f=s.field(vec,3)
visual.sphere()



#get and display points in field
for i in range(f.x):
    for j in range(f.y):
        for k in range(f.z):
            pt=f.getp(i,j,k)
            visual.sphere(pos=[pt.x,pt.y,pt.z],radius=0.5)




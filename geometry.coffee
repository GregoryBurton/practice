
class Vector
  constructor: (x=0,y=0,z=0)->
    if (x?.constructor?.name) and (x?.constructor?.name is @constructor.name)
      @constructorFromVector(x)
    else
      @constructorFromXYZ(x,y,z)
  constructorFromVector: (V)->
    @V = (v for v in V.V)
  constructorFromXYZ: (x,y,z)->
    @V = [x,y,z]
  Print: ()->
    console.log(@V)
  plus: (V)->
    returnval = new Vector()
    returnval.V = ( (@V[i]+V.V[i]) for i in [0...@V.length] )
    return returnval
  minus: (V)->
    returnval = new Vector()
    returnval.V = ( (@V[i]-V.V[i]) for i in [0...@V.length] )
    return returnval
  dot: (V)->
    returnval = 0
    (returnval += (@V[i]*V.V[i])) for i in [0...@V.length]
    return returnval
  scaled: (S)->
    returnval = new Vector()
    returnval.V = ( (@V[i]*S) for i in [0...@V.length] )
    return returnval
  cross: (V)->
    returnval = new Vector()
    returnval.V[0] = @V[1]*V.V[2] - @V[2]*V.V[1]
    returnval.V[1] = -(@V[0]*V.V[2] - @V[2]*V.V[0])
    returnval.V[2] = @V[0]*V.V[1] - @V[1]*V.V[0]
    return returnval
  norm: ()->
    returnval = Math.sqrt(@dot(this))
    return returnval
  normalized: ()->
    n = @norm()
    #console.log(n)
    if n>0
      returnval = @scaled(1.0/n)
      #console.log(returnval.V)
      return @scaled(1.0/n)
    else
      return new Vector(this)
  normalize: ()->
    n = @norm()
    if n>0
      (@V[i]/=n) for i in [0...@V.length]
  isLessThan: (V)->
    for i in [0...3]
      if @V[i]<V.V[i]
        return true
      else if @V[i]>V.V[i]
        return false
    return false
  allElementsEqual: (V)->
    for i in [0...3]
      if @V[i] isnt V.V[i]
        return false
    return true
  isEqualTo: (V)->
    return (this is V) or @allElementsEqual(V)


class Edge
  constructor: (@V0, @V1)->
    @V = @V1.minus(@V0)
    if @V1.isLessThan(@V0)
      @UndirectedV = @V0.minus(@V1)
    else
      @UndirectedV = @V
  mag: ()->
    return @V.norm()
  isNeighbor: (E)->
    return @undirectedHash()==(E.undirectedHash())
  getPortionAboveZ: (zslice)->
    if (not zslice?) or (@V0.V[2]>=zslice and @V1.V[2]>=zslice)
      return this
    if(@V0.V[2]<zslice and @V1.V[2]<zslice)
      return undefined
    if @UndirectedV is @V
      Vstart = @V0
      Vend = @V1
    else
      Vstart = @V1
      Vend = @V0
    dz = zslice - Vstart.V[2]
    VSlice = Vstart.plus(@UndirectedV.scaled(dz/@UndirectedV.V[2]))
    if @V0.V[2]<zslice
      return new Edge(VSlice, @V1)
    else
      return new Edge(@V0, VSlice)
  undirectedHash: ()->
    if @V0.id<@V1.id
      return @V0.id+','+@V1.id
    else
      return @V1.id+','+@V0.id
  consolelog: ()->
    console.log('['+@V0.V.join(',')+'] -> ['+@V1.V.join(',')+']')


class Triangle
  constructor: (V0, V1, V2)->
    @Verts = [V0, V1, V2]
    @Edges = ((new Edge(@Verts[i], @Verts[(i+1)%3])) for i in [0...3])
    @recomputeNormal = true
    @Neighbors = []
    for edge in @Edges
      edge.Triangle = this
  normal: ()->
    @N = {}
    if @recomputeNormal is true
      #console.log('got here')
      #@N = @Edges[0].V.cross(@Edges[1].V).normalized()
      @N = @Edges[0].V.cross(@Edges[1].V)
      #console.log(@N.V)
      @N = @N.normalized()
      #console.log(@N.V)
      recomputeNormal = false
    #console.log(@N.V)
    return @N
  isVisible: (zslice)->
    #console.log(@normal().V)
    if @normal().V[2]<=0.0
      return false
    #Triangle points up, so see if any vertex is above the line
    if not zslice?
      return true
    for v in @Verts
      if v.V[2]>zslice
        return true
    return false
  getVisibleBoundaryEdges: (zslice, returnval=[])->
    if (@isVisible(zslice) is false)
      return returnval
    if (@Neighbors.length isnt 3)
      console.log('Neighbors not initialized!')
      return returnval
    for edge, i in @Edges when (@Neighbors[i].isVisible(zslice) is false)
      edgeabove =  edge.getPortionAboveZ(zslice)
      if edgeabove isnt undefined
        returnval.push(edgeabove)
    if zslice?
      #Also want to slice this triangle with zslice
      anybelow = false
      for v in @Verts when v.V[2]<zslice
        anybelow = true
        break
      if anybelow
        #We do need to slice this triangle. Find a vertex above the line
        indexabove=0
        for v,i in @Verts when v.V[2]>=zslice
          indexabove = i
          break
        #Now find the first index below
        for i in [(indexabove+1)...(indexabove+3)] when @Verts[i%3].V[2]<zslice
          indexbelow = i%3
          indexabove = indexbelow-1
          if indexabove<0
            indexabove += 3
          break
        #Slice that edge
        edgeabove = @Edges[indexabove].getPortionAboveZ(zslice)
        V0 = edgeabove.V1
        #Now find the next index above
        for i in [(indexbelow+1)...(indexbelow+3)] when @Verts[i%3].V[2]>=zslice
          indexabove = i%3
          indexbelow = indexabove-1
          if indexbelow<0
            indexbelow += 3
          break
        #Slice that edge
        edgeabove = @Edges[indexbelow].getPortionAboveZ(zslice)
        V1 = edgeabove.V0
        #Create a new edge from V0 to V1
        sliceedge = new Edge(V0, V1)
        returnval.push(sliceedge)


    return returnval

class Polyhedron
  constructor: (@Verts, @Triangles)->
    #do nothing
    [scale, offset] = @getScaleAndOffset()
    #//Scale and offset the vertex
    for vert,i in @Verts
      vert.V = ( ( (500.0/scale)*(v-offset[j]) ) for v,j in vert.V)
  Print: ()->
    console.log('hello')
  getBounds: ()->
    if not(@Verts?) or (@Verts.length<1)
      return []
    #Initialize min and max to the values of the first vertex
    bounds = ([v, v] for v in @Verts[0].V)
    #Go through all of the vertices
    for vert in @Verts
      #Reset the min for all indices
      (bounds[i][0]=v) for v,i in vert.V when bounds[i][0]>v
      #Reset the max for all indices
      (bounds[i][1]=v) for v,i in vert.V when bounds[i][1]<v
    return bounds
  getScaleAndOffset: (bounds)->
    if not(bounds?)
      bounds = @getBounds()
    if bounds.length<1
      return [1, [0,0,0]]
    offset = (.5*(v[0]+v[1]) for v in bounds)
    #Max of an array: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/max
    scale = Math.max.apply(null, (Math.abs(v[0]-v[1]) for v in bounds) )
    return [scale, offset]
  getVisibleBoundaryEdges: (zslice=-1000)->
    #console.log('Getting visible boundary above '+zslice)
    returnval = []
    tri.getVisibleBoundaryEdges(zslice, returnval) for tri in @Triangles
    #e.consolelog() for e in returnval
    return returnval



###
V0 = new Vector(1,2,1)
V0.Print()
V1 = new Vector(V0)
V1.Print()
V2 = V0.plus(V1)
V2.Print()
###

###
STLParserASCII = (filetext)->
  outtext = []
  #Get rid of any \r values
  filetext = filetext.replace('\r','')
  #Split into lines and skip the first and last lines
  filetext = filetext.split('\n')[1...-1]
  for i in[0...filetext.length] by 7
    #The lines with vertex data are the 3rd through 5th inclusive
    vertlines = filetext[(i+2)..(i+4)]
    for v in vertlines
      #Grab the nonwhitespace tokens
      xyz = v.match(/\S+/g)[1..-1];
      outtext.push(xyz.join(',')+'<br>')
  return outtext.join('')
###

###
STLParserASCII = (filetext)->
  outtext = []
  #Get rid of any \r values
  filetext = filetext.replace('\r','')
  #Split into lines and skip the first and last lines
  uniqueVerts = {}
  filetext = filetext.split('\n')[1...-1]
  for i in[0...filetext.length] by 7
    #The lines with vertex data are the 3rd through 5th inclusive
    vertlines = filetext[(i+2)..(i+4)]
    for v in vertlines
      #Grab the nonwhitespace tokens
      xyz = (parseFloat(s).toFixed(10) for s in v.match(/\S+/g)[1..-1]);
      newkey = xyz.join(', ')
      if (newkey of uniqueVerts)
        uniqueVerts[newkey]++
      else
        uniqueVerts[newkey] = 1
      #outtext.push(xyz.join(', ')+'<br>')
  for key,value of uniqueVerts
    outtext.push(key+': '+value+'<br>')
  return outtext.join('')
###

STLParserASCII = (filetext)->
  outtext = []
  uniqueVerts = {}  #container of the unique vertices of this model
  numUniqueVerts = 0
  triangles=[]
  #Get rid of any \r values
  filetext = filetext.replace('\r','')
  #Split into lines and skip the first and last lines
  filetext = filetext.split('\n')[1...-1]
  for i in[0...filetext.length] by 7
    triverts = [] #array of the vertices of this triangle
    #The lines with vertex data are the 3rd through 5th inclusive
    for v in filetext[(i+2)..(i+4)]
      #Grab the nonwhitespace tokens
      xyz = (parseFloat(s) for s in v.match(/\S+/g)[1..-1]);
      #Convert the xyz into a string for the map, keep 10 decimal places to be pretty darn sure we aren't truncating
      newkey = (v.toFixed(10) for v in xyz).join(', ')
      #Pull corresponding vertex out of the map, or create a new one
      if (newkey of uniqueVerts)
        triverts.push(uniqueVerts[newkey])
      else
        newvert = new Vector(xyz[0], xyz[1], xyz[2])
        newvert.id = numUniqueVerts
        uniqueVerts[newkey] = newvert
        numUniqueVerts++
        triverts.push(newvert)
    if triverts.length is 3
      triangles.push(new Triangle(triverts[0],triverts[1],triverts[2]))
  #Define edge neighbors using the undirectedHash
  uniqueEdges = {}
  for tri in triangles
    for edge in tri.Edges
      newkey = edge.undirectedHash()
      if not (newkey of uniqueEdges)
        uniqueEdges[newkey] = []
      uniqueEdges[newkey].push(edge)
  for key,halfEdgePair of uniqueEdges
    if halfEdgePair.length isnt 2
      console.log('Non-water tight stl file!')
    halfEdgePair[0].NeighborEdge = halfEdgePair[1]
    halfEdgePair[1].NeighborEdge = halfEdgePair[0]
  for tri in triangles
    for edge,i in tri.Edges
      tri.Neighbors[i] = edge.NeighborEdge.Triangle

  alluniqueverts = (vert for key,vert of uniqueVerts)
  polyhedron = new Polyhedron(alluniqueverts, triangles)
  ###
  #Create an output
  for tri in triangles
    text = (x.id for x in tri.Verts).join(', ')
    outtext.push(text+'<br>')
  return outtext.join('')
  ###
  return polyhedron


window.STLParserASCII = STLParserASCII
console.log('geometry.coffee loaded')








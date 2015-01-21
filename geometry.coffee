
loopsFromPixelCoords = (pixels, rows, cols)->
  #Search order is Down, Left, Up, Right relative to current pixel
  search_delta = [-cols, -1, cols, 1]
  search_index = 0
  #Create a hash of the edge pixels, also find the minimum pixel
  minpel = rows*cols #Initialize to something large
  edgepixels = []
  for pel in pixels
    edgepixels[pel] = 1
    if minpel>pel
      minpel=pel
  #Start from minpel and search around
  foundpelcount = 0
  startpel = minpel
  loops = []
  while foundpelcount<pixels.length
    #Create a loop starting from startpel
    newloop = []
    pel = startpel
    while edgepixels[pel] is 1
      edgepixels[pel] = 0
      newloop.push(pel)
      foundpelcount++
      #Find the next pel
      for i in [0...4]
        j = pel + search_delta[(search_index+i)%4]
        if edgepixels[j]?
          #Found the next edge
          search_index = (search_index+3)%4  #Set the search for the next pixel
          pel = j
    loops.push(newloop)
    #Find a new startpel
    for pel,val of edgepixels when (val is 1)
      startpel = pel
      break
  return loops


window.loopsFromPixelCoords = loopsFromPixelCoords

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
  hashString: ()->
    newkey = (v.toFixed(10) for v in @V).join(', ')
    return newkey
  printString: ()->
    outstring = '(' + (v.toFixed(2) for v in @V).join(', ') + ') '
    return outstring


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

convertCollectionOfEdgesIntoLoops = (edges)->
  #This function assumes that all loops made of the collection of edges
  #are crack free. It also assumes that it is safe to return a self-intersecting loop
  verts2edges = {}
  for edge in edges
    key = edge.V0.hashString()
    if key of verts2edges
      verts2edges[key].push(edge)
    else
      verts2edges[key] = [edge]
  loops = []
  for key,keyedges of verts2edges
    for edge in keyedges when not edge.visited?
      edge.visited = true
      edgeloop = [ edge ]
      nextkey = edge.V1.hashString()
      anyappends = true #In case there is a crack and we can't ever get back to the beginning. That should never happen.
      while (nextkey isnt key) and (anyappends is true)
        anyappends = false
        nextedges = verts2edges[nextkey]
        #Just grab the first unvisited edge, when there is more than one it doesn't matter which you pick
        for nextedge in nextedges when not nextedge.visited?
          nextedge.visited = true
          edgeloop.push(nextedge)
          nextkey = nextedge.V1.hashString()
          anyappends = true
          break
      loops.push(edgeloop)
  for edge in edges
    delete edge.visited
  return loops

printEdgeLoop = (edgeloop)->
  console.log('0: '+edgeloop[0].V0.printString())
  for e,i in edgeloop
    console.log(i+': '+edgeloop[i].V1.printString())
  return

connectColinearSegmentsInEdgeLoop = (edgeloop)->
  if edgeloop.length<2
    return edgeloop
#  console.log(edgeloop)
  newloop = [edgeloop[0]]
  lastedge = edgeloop[0]
  for i in [1...edgeloop.length]
    testedge = new Edge(lastedge.V0, edgeloop[i].V1)
    #Check if this edge is colinear with lastedge
    if( lastedge.V.normalized().dot( testedge.V.normalized() ) > (1.0-1e-6) )
      lastedge = testedge
      newloop[newloop.length-1] = testedge
    else
      newloop.push(edgeloop[i])
      lastedge = edgeloop[i]
  #At this point, it is possible that the last edge and first edge are colinear
  testedge = new Edge(lastedge.V0, newloop[0].V1)
  if( lastedge.V.normalized().dot( testedge.V.normalized() ) > (1.0-1e-6) )
#    console.log('combining last and first')
    #Change the start point of newloop[0] and pop the last edge
    newloop[0] = testedge
    newloop.pop()
#  console.log('length was ' + edgeloop.length + ' now is ' + newloop.length)
#  printEdgeLoop(edgeloop)
#  printEdgeLoop(newloop)
  return newloop


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
  sliceAtZ: (zslice)->
    #Return the edge representing the zplane cut through this triangle
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
      return sliceedge
    else
      return null
  getVisibleBoundaryEdges: (zslice, returnval=[])->
    if @visited?
      #We've already visited this triangle
      return returnval
    #We're visiting this triangle
    @visited = true
    if (@isVisible(zslice) is false)
      return returnval
    if (@Neighbors.length isnt 3)
      console.log('Neighbors not initialized!')
      return returnval
    if zslice?
      #Want to slice this triangle with zslice
      sliceedge = @sliceAtZ(zslice)
      if sliceedge
        returnval.push(sliceedge)
    #Get edges that bound visible/nonvisible transitions
    for edge, i in @Edges
      if (@Neighbors[i].isVisible(zslice) is false)
        #This edge is the transition from visible to nonvisible
        #Get the portion of the edge that is above the zplane
        edgeabove =  edge.getPortionAboveZ(zslice)
        if edgeabove isnt undefined
          returnval.push(edgeabove)
      else
        #The neighbor is visible, so just recurse into that triangle to find more boundaries
        @Neighbors[i].getVisibleBoundaryEdges(zslice, returnval)
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
  verifyNoCrackInLoop: (returnval)->
    #Verify that the edge loops have no cracks  Vert.hashString()
    vertexHashes = {}
    for edge in returnval
      edgeverts = [edge.V0, edge.V1]
      for vert in edgeverts
        newkey = vert.hashString()
        if (newkey of vertexHashes)
          #Increase counts for that vertex
          vertexHashes[newkey][1] += 1
          #Should replace the outbound edge vertex with the stored vertex
        else
          #Store the vertex with this hash and initialize the count to 1
          vertexHashes[newkey] = [vert, 1]
    #The way we verify that the edge loops have no cracks is to check that the
    #count for each vertex is even. We know this because in a complete loop there
    #should be an edge in and an edge out from each vertex. Given that you can
    #many triangles fanning out from a single vertex, the number of edges in/out of that
    #vertex may be greater than 2, but it must be an even number.
    thereIsACrack = false
    for key,value of vertexHashes
      console.log(key+': '+value[1])
      if value[1]%2 isnt 0
        console.log('CRACK IN LOOP!')
        thereIsACrack = true
    if thereIsACrack
      console.log('There is a crack!')
    else
      console.log('No crack in any loop!')
    return not(thereIsACrack)
  getVisibleBoundaryEdges: (zslice=-1000)->
    #console.log('Getting visible boundary above '+zslice)
    returnval = [] #a collection of collections of edges (each collection comes from a connected set of triangles)
    for tri,i in @Triangles when (tri.visited? is false)
      collection = []
      tri.getVisibleBoundaryEdges(zslice, collection)
      if collection.length>0
        #Grab all of the loops as separately
        (returnval.push(edgeloop) for edgeloop in convertCollectionOfEdgesIntoLoops(collection) )
        #console.log('collection '+i+' = '+collection.length)
      tri.visited = true
    #e.consolelog() for e in returnval
    #@verifyNoCrackInLoop(returnval)
    for tri in @Triangles
      delete tri.visited
    #console.log('Num collections = '+returnval.length)
    return returnval
  getVisibleBoundaryEdgesPlanarized: (zslice=0, offsetamount=20)->
    boundaryloops = @getVisibleBoundaryEdges(zslice)
    outputloops = []  #This will be some subset of the edges of the boundaryloops, flattened to zslice
    uniqueVerts = {}
    for edgeloop in boundaryloops when (edgeloop.length>2)
      flatloop = []
      for edge in edgeloop
        V = [new Vector(edge.V0.V[0],edge.V0.V[1],zslice), new Vector(edge.V1.V[0],edge.V1.V[1],zslice)];
        for v,i in V
          newkey = v.hashString()
          if v.newkey of uniqueVerts
            V[i] = uniqueVerts[newkey]
          else
            uniqueVerts[newkey] = v
        newedge = new Edge(V[0], V[1])
        if newedge.mag()>0
          flatloop.push(newedge)
#      cleanedloops = connectColinearSegmentsInEdgeLoop(flatloop)
#      outputloops.push( cleanedloops )
      outputloops.push( flatloop )
#    return outputloops
    offsetloops = clipperOffsetLoops(outputloops, offsetamount, zslice)
    return offsetloops

vertex2ClipperVertex = (vert, scale=1.0)->
  cvert = {}
  cvert.X = vert.V[0]*scale
  cvert.Y = vert.V[1]*scale
  return cvert

clipperOffsetLoops = (loops, offsetamount, z)->
  #Create clipper loops out of each loop
  clipperscale = 1e5
  cloops = []
  for edgeloop in loops
#    console.log(edgeloop)
    newloop = [ vertex2ClipperVertex(edgeloop[0].V0, clipperscale) ]
    for edge in edgeloop
      newloop.push(vertex2ClipperVertex(edge.V1, clipperscale))
    cloops.push(newloop)
  ###
  #First, union all of the loops
  cpr = new ClipperLib.Clipper()
  cpr.AddPaths(cloops, ClipperLib.PolyType.ptSubject, true) #true means "closed paths" instead of open paths
  solution_paths = new ClipperLib.Paths()
  clipType = ClipperLib.ClipType.ctUnion
  subject_fillType = ClipperLib.PolyFillType.pftNonZero
  clip_fillType = ClipperLib.PolyFillType.pftNonZero
  succeeded = cpr.Execute(clipType, solution_paths, subject_fillType, clip_fillType)
  ###
  solution_paths = ClipperLib.Clipper.SimplifyPolygons(cloops, ClipperLib.PolyFillType.pftNonZero)
  cleandelta = 1e-3
  solution_paths = ClipperLib.Clipper.CleanPolygons(solution_paths, cleandelta * clipperscale)

  #Create clipper offset object
  miterLimit = 2
  arcTolerance = 1e-3*clipperscale #0.25
  co = new ClipperLib.ClipperOffset(miterLimit, arcTolerance)
  co.AddPaths(solution_paths, ClipperLib.JoinType.jtRound, ClipperLib.EndType.etClosedPolygon)
  offsetted_polytree = new ClipperLib.PolyTree()
#  co.Execute(offsetted_polytree, offsetamount*clipperscale)
#  console.log(offsetted_polytree)
#  window.polytree = offsetted_polytree
  offsetted_paths = new ClipperLib.Paths()
  co.Execute(offsetted_paths, offsetamount*clipperscale)
  #console.log(offsetted_paths)
  #Create my loops out of clipper loops
  scale_down = 1e-5
  loopsout = []
  for cloop in offsetted_paths
    loopout = []
    l = cloop.length
    v0 = new Vector(cloop[0].X*scale_down, cloop[0].Y*scale_down, z)
    for vert,i in cloop
      if i is 0
        continue
      v1 = new Vector(vert.X*scale_down, vert.Y*scale_down, z)
      newedge = new Edge(v0, v1)
      loopout.push(newedge)
      v0 = v1
    if loopout[0].V0 isnt loopout[loopout.length-1].V1
      newedge = new Edge(loopout[loopout.length-1].V1, loopout[0].V0)
      loopout.push(newedge)
    loopsout.push(loopout)
  return loopsout

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
window.DistPel = DistPel

###
loop 1
embedded:11 (250.00, 171.70, 0.00)
2embedded:11 (-250.00, 171.70, 0.00)
2embedded:11 (-250.00, -171.70, 0.00)
2embedded:11 (250.00, -171.70, 0.00)
embedded:11 (250.00, 171.70, 0.00)
embedded:14 loop 2
embedded:11 (-113.13, -100.97, 0.00)
2embedded:11 (-113.13, 85.66, 0.00)
2embedded:11 (162.18, 85.66, 0.00)
2embedded:11 (162.18, -100.97, 0.00)
embedded:11 (-113.13, -100.97, 0.00)
embedded:14 loop 3
embedded:11 (76.14, 17.23, 0.00)
2embedded:11 (-62.89, 17.23, 0.00)
2embedded:11 (-62.89, -63.29, 0.00)
2embedded:11 (76.14, -63.29, 0.00)
embedded:11 (76.14, 17.23, 0.00)
###





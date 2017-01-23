#!/usr/bin/ruby

require 'matrix'

require 'cmath'
include CMath

require 'daru'
require 'gnuplotrb'
include Daru
include GnuplotRB

def dist_uh(x, y)
  Math::acosh(1.0 + (x-y).abs**2.to_f / (2 * x.imag * y.imag))
end

def dist_pd(x, y)
  dist_uh(to_uh(x), to_uh(y))
end

def to_uh(x)
  c_uh = Matrix[[1, 1i], [1i, 1]]
  (c_uh * Vector[x, 1]).inject(:/)
end

def to_pd(x)
  c_pd = Matrix[[1, 1i], [1i, 1]].inverse
  (c_pd * Vector[x,1]).inject(:/)
end

def dirichlet_domain(p, generators, depth = 3)
  orbit = get_orbit(p, generators, depth)

  circles = []

  plot = get_boundary_plot

  orbit.each do |q|
    pt = DataFrame.new({x: [q.real], y: [q.imag]})
    plot = (plot << [pt, using: '2:3'])

    circ = get_bisector(p, q)

    if circ.r > 1.0e5
      puts "Warning: Circle too large"
      puts q
      next
    end
    
    admissible = true
    lc = []
    circles.each do |c|
      i = triple_inclusion(p, c, circ)
      if i == 0
        lc << c
      elsif i == 1
        lc << c
        admissible = false
      end
    end

    lc << circ if admissible
    circles = lc
  end
  
  circles.each do |circ|
    df = generalized_circle_points(circ)
    l = std_circle_plot(df)
    plot = (plot << l)
  end
  pt = DataFrame.new({x: [p.real], y: [p.imag]})
  plot = (plot << [pt, using: '2:3'])
  plot
end

def plot_tesselation(tess, pts = true)
  plot = get_boundary_plot

  all_circles = []
  
  tess.each do |p, circles|

    plot = plot << plot_pt(p) if pts
    circles -= all_circles
    circles.each do |circ|
      next unless is_geodesic?(circ)
      plot = plot << plot_circle(circ)
    end
    all_circles += circles
  end
  plot
end

def get_tesselation(p, generators, depth = 3, orbit = nil)
  orbit, root = get_orbit(p, generators, depth) if orbit == nil

  tess = {}

  tess[p] = get_intersection(p, orbit)

  stack = root.children
  
  until stack.empty? do
    child = stack.pop
    g = child.lastg
    stack += child.children
    tess[child.p] = []

    tess[p].each do |circ|
      tess[child.p] << translate_by(circ, g)
    end
  end

  tess
end

def get_intersection(p, orbit)
  circles = []

  orbit.each do |q|
    circ = get_bisector(p, q)

    if not circ.r.nil? and circ.r > 1.0e5
      puts "Warning: Circle too large"
      puts q
      next
    end

    admissible = true
    lc = []
    circles.each do |c|
      i = triple_inclusion(p, c, circ)
      if i == 0
        lc << c
      elsif i == 1
        lc << c
        admissible = false
      end
    end

    lc << circ if admissible
    circles = lc
  end
  circles
end

def plot_pt(p)
  pt = DataFrame.new({x: [p.real], y: [p.imag]})
  [pt, using: '2:3']
end

def plot_circle(c)
  df = generalized_circle_points(c)
  std_circle_plot(df)
end

Caley = Struct.new(:p, :lastg, :length, :children)

def get_orbit(p, generators, depth, eps = 0.001)

  all = generators.map do |g|
    g.inverse
  end

  all += generators

  root = Caley.new(p, Matrix[[1,0], [0,1]], 0, [])

  stack = [root]

  orbit = []

  return orbit if depth == 0

  until stack.empty? do
    vertex = stack.pop

    all.each do |g|
      next if g.inverse == vertex.lastg

      q = g * Vector[vertex.p, 1]
      q = q[0] / q[1]

      skip = false
      orbit.each do |o|
        skip = true if (q - o).abs < eps
      end
      next if skip or (q - p).abs < eps

      orbit << q

      unless vertex.length + 1 >= depth
        child = Caley.new(q, g, vertex.length + 1, [])
        stack << child
        vertex.children << child
      end
    end
  end

  [orbit, root]
end

Circle = Struct.new(:r, :c)

# in this case we know that if there is a line, it will go through the origin.
def get_bisector(p, q)

  qa = q.abs
  pa = p.abs
  
  a = qa**2 - pa**2
  b = (1 - pa**2) * q.conjugate - (1 - qa**2) * p.conjugate
  c = b.conjugate
  d = a

  if a == 0
    return Circle.new(nil, 1.0i*(p-q))
  end

  r = Math::sqrt((b * c).real / a**2 - d / a)
  c = c / a
  
  Circle.new(r, c)
end

def generalized_circle_points(circ, step = 0.01)

  # line
  if circ.r.nil?

    v = circ.c / circ.c.abs

    x = []
    y = []
    
    (-1..1).step(step).each do |t|
      
      x << t * v.real
      y << t * v.imag
    end
    
  # circle
  else
    x = []
    y = []

    alpha = Math::acos(circ.r / circ.c.abs)

    gamma = circ.c.arg + PI
    
    (-alpha..alpha).step(step).each do |t|

      p = circ.c + circ.r * exp((t + gamma) * 1i)

      x << p.real
      y << p.imag
      
    end
  end

  DataFrame.new({x: x, y: y})
end

def get_boundary(step = 0.01)

  x = []
  y = []

  (0..2*PI).step(step).each do |t|

    p = exp(t * 1i)

    x << p.real
    y << p.imag
  end

  DataFrame.new({x: x, y: y})
end

def compute_schottky_element(c1, c2)

  p1 = get_triple(c1)
  p2 = get_triple(c2).reverse

  m1 = cross_ratio(p1)
  m2 = cross_ratio(p2).inverse

  m2 * m1  
end

def get_triple(circ)
  if circ.r.nil?
    v = circ.c / circ.c.abs
    p1 = [v, 0, -v]
  else
    alpha = Math::acos(circ.r / circ.c.abs)
    gamma = circ.c.arg + PI
    p1 = [circ.c + circ.r * exp((alpha + gamma) * 1i)]
    p1 << circ.c + circ.r * exp(gamma * 1i)
    p1 << circ.c + circ.r * exp((-alpha + gamma) * 1i)
  end

  p1
end

def cross_ratio(zs)

  f = (zs[1] - zs[2]) / (zs[1] - zs[0])

  Matrix[[f, -f*zs[0]], [1, -zs[2]]]
end

def apply_schottky(df, m)

  x = []
  y = []

  df.each_row do |r|
    z = Complex(r[:x], r[:y])
    v = (m * Vector[z, 1]).inject(:/)
    x << v.real
    y << v.imag
  end

  DataFrame.new({x: x, y: y})
end

def get_boundary_plot
  bounding_circle = get_boundary
  Plot.new([bounding_circle, with: 'lines', using: '2:3', lt: 'rgb "black"'],
           xrange: -1..1,
           yrange: -1..1,
           border: '0',
           ytics: nil,
           xtics: nil
          )
end

def std_circle_plot(df)
  [df, with: 'lines', using: '2:3', lt: 'rgb "blue"']
end

def find_neighbours(p, orbit, len)

  dist = orbit.map do |o|
    dist_pd(p, o)
  end
  
  df = DataFrame.new({points: orbit, dist: dist})
  min = df[:dist].min
  index = df[:dist].map do |d|
    d < 2*min
  end
  df[:points].where(index).to_a
end

# Handle lines
def in_circ?(p, c)
  return false if c.r.nil?
  (p - c.c).abs < c.r
end

def circ_in_circ?(c1, c2, eps = 0.01)
  return false if c1.r.nil? or c2.r.nil?
  (c1.c - c2.c).abs + c1.r < c2.r + eps
end

def triple_inclusion(p, c1, c2)

  # TODO: handle lines
  if c1.r.nil? or c2.r.nil?
    return 0
  end
  
  dist = (c1.c - c2.c).abs
  r_max = [c1.r, c2.r].max
  r_min = [c1.r, c2.r].min

  if dist + r_min >= r_max
    return 0
  else
    if in_circ?(p, c1) and in_circ?(p, c2)
      return c1.r < c2.r ? 1 : 2
    elsif (not in_circ?(p, c1)) and (not in_circ?(p, c2))
      return c1.r > c2.r ? 1 : 2
    else
      return 0
    end
  end
end

def translate_by(circ, g)
  if circ.r.nil?
    a = 0.0
    b = circ.c*1.0i
    c = b.conjugate
    d = 0.0
  else
    a = 1.0
    b = - circ.c.conjugate
    c = b.conjugate
    d = circ.c.abs**2 - circ.r**2
  end

  h = proj_inv(g)

  at = (a*h[0,0]**2 + (b+c)*h[0,0]*h[1,0] + d*h[1,0]**2).real
  bt = a*h[0,0]*h[0,1] + b*h[0,0]*h[1,1] + c*h[0,1]*h[1,0] + d*h[1,0]*h[1,1]
  ct = bt.conjugate
  dt = (a*h[0,1]**2 + (b + c)*h[0,1]*h[1,1] + d*h[1,1]**2).real

  if at == 0.0
    puts "Warning! Nonzero interesct in translated circle." if d != 0.0
    return Circle.new(nil, 1.0i*bt)
  else
    return Circle.new(Math::sqrt((bt*ct).real/at**2 - dt/at), - ct/at)
  end
end

def proj_inv(m)
  Matrix[[m[1,1], -m[0,1]], [-m[1,0], m[0,0]]]
end

def is_geodesic?(c, eps = 0.01)
  (c.c.abs**2 - 1.0 - c.r**2)**2 < eps**2
end

c1 = Circle.new(1.0, Complex(Math::sqrt(2), 0.0))
c2 = Circle.new(1.0, Complex(0.0, Math::sqrt(2)))
c3 = Circle.new(0.5, Complex(-Math::sqrt(1.25), 0.0))
c4 = Circle.new(0.5, Complex(0.0, -Math::sqrt(1.25)))

a = compute_schottky_element(c1, c3)
b = compute_schottky_element(c2, c4)

gen = [a,b]
p = 0.0

puts "load 'fd-plot.rb'; plot = dirichlet_domain(p, gen, 2); plot.to_png('test.png'," +
     "size: [700, 700]); nil"

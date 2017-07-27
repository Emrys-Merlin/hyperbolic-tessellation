# frozen_string_literal: true
require_relative 'matrix.rb'
require_relative 'circle.rb'

# This class contains part of the Caley graph of a group of moebius transformations.
class Graph
  include Enumerable
  attr_reader :depth, :root, :gen

  def initialize(p, gen, depth)
    @prng = Random.new
    @root = Node.new(p, Matrix[[1, 0], [0, 1]], Matrix[[1, 0], [0, 1]],
                     0, [], [], {})
    @depth = 0
    @gen = gen.map(&:proj_inv)
    @gen += gen
    deepen(depth)
  end

  def deepen(depth, eps = 0.001)
    stack = [@root]

    until stack.empty?
      node = stack.pop

      next if node.depth >= depth

      if node.depth >= @depth
        @gen.each do |g|
          next if g.proj_inv == node.prefix

          q = g.moebius(node.p)

          skip = false
          each do |n|
            o = n.p
            if (q - o).abs < eps
              skip = true
              break
            end
          end

          node.children << Node.new(q, g * node.g, g, node.depth + 1,
                                    [], [], {})
        end
      end

      stack += node.children
    end

    @depth = depth
  end

  # Computes the Dirichlet fundamental domain around the point center
  def get_intersection(center = @root, nr_tests = 500)
    circles = []

    each do |node|
      next if node == center

      circles << Circle.get_bisector(center.p, node.p)

      # admissible = true
      # lc = []
      # circles.each do |c|
      #   i = c.triple_inclusion(center.p, circ)
      #   if i.zero?
      #     lc << c
      #   elsif i == 1
      #     lc << c
      #     admissible = false
      #     # break
      #   end
      # end
      # lc << circ if admissible
      # circles = lc
    end

    circles.uniq!
    admissibles = [false] * circles.length
    get_ud_pts(nr_tests).each do |pt|
      counter = 0
      index = -1
      circles.length.times do |i|
        if circles[i].separates(center.p, pt)
          counter += 1
          index = i
        end
      end
      admissibles[index] = true if counter == 1
    end

    center.border = []

    admissibles.length.times do |i|
      center.border << circles[i] if admissibles[i]
    end

    # center.border = circles
    center.angles = get_angles(center)
  end

  # Translate the Dirichlet fundamental domain around center to all
  # its orbit points (up to depth)
  def get_tessellation(center = @root)
    get_intersection(center)

    each do |node|
      next if node == center

      node.border = center.border.map do |circ|
        circ.translate_by(node.g)
      end

      node.angles = get_angles(node)
    end
  end

  # a pre-order depth first search implementation to visit all nodes
  def each
    queue = [@root]
    until queue.empty?
      node = queue.shift
      queue += node.children
      yield node
    end
  end

  def nodes
    nodes = []
    each do |node|
      nodes << node
    end
    nodes
  end

  def get_angles(node)
    circles = node.border
    angles = Hash.new { |hsh, key| hsh[key] = [] }
    # Find intersections between each pair of circles c1, c2
    (0...circles.length).each do |i|
      c1 = circles[i]
      ((i + 1)...circles.length).each do |j|
        c2 = circles[j]
        pts = c1.get_intersection(c2)

        next if pts.nil?
        # Get the angle for each point lying in the unit disk
        pts.each do |pt|
          next if pt.abs >= 1.0
          angles[c2] << (pt - c2.c).angle.round(6)
          angles[c2].uniq!
          angles[c1] << (pt - c1.c).angle.round(6)
          angles[c1].uniq!
        end
      end

      # If c1 does not have any (so far: enough) intersections, get the intersections with
      # the unit disk to find the right segment.
      next if angles[c1].length >= 2
      angles[c1] = [0, 2 * Math::PI]

      pts = c1.get_intersection(Circle.new(c: 0.0+0.0i, r: 1.0))
      next if pts.nil?
      angles[c1] = []
      pts.each do |pt|
        angles[c1] << (pt - c1.c).angle
      end

      # angles[c1].map! { |a| a.round(8) }.uniq!
    end
    angles
  end

  def get_ud_pts(n)
    result = []

    n.times do
      loop do 
        x = @prng.rand(-1.0..1.0)
        y = @prng.rand(-1.0..1.0)
        z = Complex(x, y)
        if z.abs < 1.0
          result << z
          break
        end
      end
    end

    result
  end
end

Node = Struct.new(:p, :g, :prefix, :depth, :children, :border, :angles)

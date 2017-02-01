# frozen_string_literal: true
require_relative 'matrix.rb'
require_relative 'circle.rb'

# This class contains part of the Caley graph of a group of moebius transformations.
class Graph
  include Enumerable
  attr_reader :depth, :root, :gen

  def initialize(p, gen, depth)
    @root = Node.new(p, Matrix[[1, 0], [0, 1]], Matrix[[1, 0], [0, 1]],
                     0, [], [])
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

          node.children << Node.new(q, g * node.g, g, node.depth + 1, [], [])
        end
      end

      stack += node.children
    end

    @depth = depth
  end

  # Computes the Dirichlet fundamental domain around the point center
  def get_intersection(center = @root)
    circles = []

    each do |node|
      next if node == center

      circ = Circle.get_bisector(center.p, node.p)

      admissible = true
      lc = []
      circles.each do |c|
        i = c.triple_inclusion(center.p, circ)
        if i.zero?
          lc << c
        elsif i == 1
          lc << c
          admissible = false
          # break
        end
      end
      lc << circ if admissible
      circles = lc
    end
    center.border = circles
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
end

Node = Struct.new(:p, :g, :prefix, :depth, :children, :border)

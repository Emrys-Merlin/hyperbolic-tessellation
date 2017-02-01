# frozen_string_literal: true
require_relative 'circle.rb'

# This module generates (symmetric or regular) schottky groups with
# a given number of generators.
module Schottky
  def self.generate_circles(gen = 1)
    alpha = Math::PI / (2 * gen)
    scale = 0.9
    r = Math.tan(alpha * scale)
    r_c = 1.0 / Math.cos(alpha * scale)

    circles = []

    gen.times do |j|
      circles << Circle.new(c: r_c * exp(2i * alpha * j), r: r)
      circles << Circle.new(c: r_c * exp(2i * alpha * j + Math::PI * 1i), r: r)
    end
    circles
  end

  def self.generate_transformations(gen = 1)
    circles = generate_circles(gen)
    generators = []
    gen.times do |j|
      first = 2 * j
      last = 2 * j + 1
      generators << circles[first].get_transformation(circles[last])
    end
    generators
  end
end

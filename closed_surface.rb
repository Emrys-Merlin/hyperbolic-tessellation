# frozen_string_literal: true
require_relative 'circle.rb'

# Generate the bounding circles and transformations for closed surfaces
# of higher genus
module ClosedSurface
  def self.generate_circles(genus = 2)
    alpha = Math::PI / (4 * genus)
    s = Math.sin(alpha)
    t = Math.cos(alpha)

    dist = 1.0 / Math.sqrt(1 - (s / t)**2)
    r = 1.0 / Math.sqrt((t / s)**2 - 1)

    circles = []
    (4 * genus).times do |j|
      circles << Circle.new(c: dist * exp(alpha * 2i * j), r: r)
    end
    circles
  end

  def self.generate_transformations(genus = 2)
    circles = generate_circles(genus)

    gen = []
    genus.times do |g|
      2.times do |j|
        first = 4 * g + j
        second = 4 * g + j + 2
        gen << circles[first].get_transformation(circles[second])
      end
    end
    gen
  end
end

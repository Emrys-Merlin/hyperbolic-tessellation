# frozen_string_literal: true
require 'complex'

# Just the standard complex class with an extra method
class Complex
  def plot_pt
    coord = [[real], [imag]]
    [coord]
  end
end

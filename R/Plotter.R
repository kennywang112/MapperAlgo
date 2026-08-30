#' Plot Mapper Result
#'
#' Visualizes the Mapper output using networkD3.
#'
#' @param Mapper Mapper object.
#' @param label Label of the data.
#' @param avg Whether coloring the nodes by average label or majority label.
#' @param use_embedding Whether to use original data for coloring (TRUE or FALSE).
#' @param legend_name Optional custom title shown on the color legend / colorbar,
#'   describing what \code{label} represents (e.g. \code{"Eigen centrality"}).
#'   Defaults to \code{"Label"}, \code{"Avg(label)"}, or \code{"Majority label"}
#'   depending on the coloring mode when left \code{NULL}.
#' @return Plot of the Mapper.
#' @importFrom networkD3 forceNetwork
#' @importFrom htmlwidgets JS
#' @importFrom rlang .data
#' @export

MapperPlotter <- function(
    Mapper, label, avg=FALSE, use_embedding=FALSE, legend_name=NULL
) {

  piv <- Mapper$points_in_vertex
  num_vertices <- Mapper$num_vertices
  vertex.size <- sapply(piv, length)

  adj_indices <- which(Mapper$adjacency == 1, arr.ind = TRUE)
  adj_indices <- adj_indices[adj_indices[, 1] < adj_indices[, 2], , drop = FALSE]

  if (nrow(adj_indices) > 0) {
    edge_weights <- apply(adj_indices, 1, function(idx) {
      length(intersect(piv[[idx[1]]], piv[[idx[2]]]))
    })
  } else {
    edge_weights <- numeric(0)
  }

  legend <- TRUE
  color_title <- "Label"

  if (use_embedding) {
    if (!is.atomic(label) || !is.null(dim(label))) {
      stop(
        "use_embedding=TRUE requires `label` to be a plain vector (one value ",
        "per Mapper vertex) -- got a ", paste(class(label), collapse = "/"),
        ifelse(is.null(dim(label)), "", paste0(" with dim ", paste(dim(label), collapse = "x"))),
        ". Did you mean to pass a single column, e.g. `embedding[, 1]`?"
      )
    }
    if (length(label) != num_vertices) {
      stop(sprintf(
        paste(
          "use_embedding=TRUE requires length(label) to equal Mapper$num_vertices.",
          "Got length(label) = %d but Mapper$num_vertices = %d.",
          "`label` must already be one value per Mapper vertex, in vertex order",
          "(e.g. the output of eigen_centrality() on igraph::graph_from_adjacency_matrix(Mapper$adjacency)) --",
          "not one value per original data point."
        ),
        length(label), num_vertices
      ))
    }
    Group_col <- label

  } else if (avg) {
    legend <- FALSE
    avg_label <- vapply(piv, function(idx) mean(label[idx], na.rm = TRUE), numeric(1))
    Group_col <- avg_label
    color_title <- "Avg(label)"

  } else {
    lab_chr <- as.character(label)
    majority <- character(num_vertices)
    for (i in seq_len(num_vertices)) {
      pts <- piv[[i]]
      if (length(pts) > 0) {
        ux <- unique(lab_chr[pts])
        majority[i] <- ux[which.max(tabulate(match(lab_chr[pts], ux)))]
      } else {
        majority[i] <- "NA"
      }
    }
    Group_col <- factor(majority)
    color_title <- "Majority label"
  }

  if (!is.null(legend_name)) {
    color_title <- legend_name
  }

  MapperNodes <- data.frame(
    Nodename = 1:num_vertices,
    Nodesize = vertex.size * 5,
    Group = Group_col
  )

  if (nrow(adj_indices) > 0) {
    MapperLinks <- data.frame(
      Linksource = adj_indices[, 1] - 1,
      Linktarget = adj_indices[, 2] - 1,
      Linkvalue = edge_weights
    )
  } else {
    MapperLinks <- data.frame(Linksource=numeric(0), Linktarget=numeric(0), Linkvalue=numeric(0))
  }

  if (is.numeric(MapperNodes$Group)) {
    rng <- range(MapperNodes$Group, na.rm = TRUE)
    colourScale <- htmlwidgets::JS(sprintf(
      "d3.scaleSequential(d3.interpolateViridis).domain([%f, %f])",
      rng[1], rng[2]
    ))
    is_continuous <- TRUE
  } else {
    colourScale <- htmlwidgets::JS("d3.scaleOrdinal(d3.schemeCategory10)")
    is_continuous <- FALSE
  }

  # forceNetwork's built-in legend is only meaningful for a categorical
  # (ordinal) colourScale -- it lists color.domain(), which for our
  # continuous d3.scaleSequential is just the 2-value [min, max]
  # interpolation range, rendered as two misleading swatches (e.g. "0"/"1").
  # So the native legend is only shown for the discrete/categorical case;
  # the continuous case gets our own colorbar below instead.
  show_native_legend <- legend && !is_continuous

  p <- forceNetwork(
    Nodes = MapperNodes,
    Links = MapperLinks,
    Source = "Linksource",
    Target = "Linktarget",
    Value  = "Linkvalue",
    NodeID = "Nodename",
    Nodesize = "Nodesize",
    Group = "Group",
    opacity = 1,
    zoom = TRUE,
    radiusCalculation = JS("Math.sqrt(d.nodesize)"),
    colourScale = colourScale,
    linkDistance = JS("function(d){ return 150 / Math.sqrt(d.value + 1); }"),
    charge = JS("function(d){ return - (60 + 2*Math.sqrt(d.nodesize)); }"),
    legend = show_native_legend
  )

  if (is_continuous) {
    pal <- viridisLite::viridis(100)
    pal_json <- jsonlite::toJSON(pal, auto_unbox = TRUE)

    p <- htmlwidgets::onRender(p, htmlwidgets::JS(sprintf(
      "function(el, x) {
           var colors = %s;
           var minv = %f, maxv = %f;
           var root = d3.select(el);
           var container = root.append('div')
             .attr('class','rd3-colorbar')
             .style('position','absolute')
             .style('right','10px')
             .style('top','10px')
             .style('padding','6px')
             .style('background','rgba(255,255,255,0.95)')
             .style('border','1px solid #ddd')
             .style('border-radius','3px')
             .style('font-family','sans-serif')
             .style('font-size','11px')
             .style('pointer-events','none');

           container.append('div').text('%s').style('margin-bottom','4px').style('font-weight','500');

           var grad = container.append('div')
             .style('width','140px')
             .style('height','12px')
             .style('border','1px solid #ccc')
             .style('background', 'linear-gradient(to right,' + colors.join(',') + ')');

           var labels = container.append('div').style('display','flex').style('justify-content','space-between').style('margin-top','4px');
           labels.append('div').text(minv.toFixed(2));
           labels.append('div').text(maxv.toFixed(2));
         }",
      pal_json, rng[1], rng[2], color_title
    )))
  } else if (show_native_legend) {
    # Give the native top-left legend (built by forceNetwork.js as one
    # <g class="legend"> per category, at a fixed translate(18, i*22+4))
    # a title, by shifting all of its groups down and adding a heading
    # above them -- one single, correctly-labelled legend instead of a
    # second disconnected box.
    p <- htmlwidgets::onRender(p, htmlwidgets::JS(sprintf(
      "function(el, x) {
           var svg = d3.select(el).select('svg');
           var legendRectSize = 18, legendSpacing = 4, titleOffset = 20;

           svg.selectAll('.legend').attr('transform', function(d, i) {
             var height = legendRectSize + legendSpacing;
             var vert = i * height + 4 + titleOffset;
             var horz = legendRectSize;
             return 'translate(' + horz + ',' + vert + ')';
           });

           svg.append('text')
             .attr('x', legendRectSize)
             .attr('y', 14)
             .attr('font-family', 'sans-serif')
             .attr('font-size', '12px')
             .attr('font-weight', '600')
             .text('%s');
         }",
      color_title
    )))
  }

  return(p)
}

#' Plot Mapper Result in Interactive 3D
#'
#' @param Mapper Mapper object.
#' @param label Label of the data.
#' @param avg Whether coloring the nodes by average label or majority label.
#' @param use_embedding Whether to use original data for coloring (TRUE or FALSE).
#' @param seed Random seed for the stochastic force-directed layout, so the plot is reproducible across calls.
#' @param legend_name Optional custom title shown on the color legend / colorbar,
#'   describing what \code{label} represents (e.g. \code{"Eigen centrality"}).
#'   Defaults to \code{"Label"}, \code{"Avg(label)"}, or \code{"Majority label"}
#'   depending on the coloring mode when left \code{NULL}.
#' @return An interactive plotly 3D plot of the Mapper graph.
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr as_edgelist
#' @importFrom plotly plot_ly add_trace layout
#' @export
MapperPlotter3D <- function(
    Mapper, label, avg=FALSE, use_embedding=FALSE, seed=42, legend_name=NULL
) {

  piv <- Mapper$points_in_vertex
  num_vertices <- Mapper$num_vertices
  vertex.size <- sapply(piv, length)

  # Build the graph once and read edges back from igraph itself, so the edge
  # order used to weight the layout and the edge order used to draw the
  # connecting lines always refer to the same edges.
  g <- igraph::graph_from_adjacency_matrix(Mapper$adjacency, mode = "undirected")
  edge_list <- igraph::as_edgelist(g, names = FALSE)

  if (nrow(edge_list) > 0) {
    edge_weights <- apply(edge_list, 1, function(idx) {
      length(intersect(piv[[idx[1]]], piv[[idx[2]]]))
    })
  } else {
    edge_weights <- numeric(0)
  }

  legend <- TRUE
  color_title <- "Label"

  if (use_embedding) {
    if (!is.atomic(label) || !is.null(dim(label))) {
      stop(
        "use_embedding=TRUE requires `label` to be a plain vector (one value ",
        "per Mapper vertex) -- got a ", paste(class(label), collapse = "/"),
        ifelse(is.null(dim(label)), "", paste0(" with dim ", paste(dim(label), collapse = "x"))),
        ". Did you mean to pass a single column, e.g. `embedding[, 1]`?"
      )
    }
    if (length(label) != num_vertices) {
      stop(sprintf(
        paste(
          "use_embedding=TRUE requires length(label) to equal Mapper$num_vertices.",
          "Got length(label) = %d but Mapper$num_vertices = %d.",
          "`label` must already be one value per Mapper vertex, in vertex order",
          "(e.g. the output of eigen_centrality() on igraph::graph_from_adjacency_matrix(Mapper$adjacency)) --",
          "not one value per original data point."
        ),
        length(label), num_vertices
      ))
    }
    Group_col <- label

  } else if (avg) {
    legend <- FALSE
    avg_label <- vapply(piv, function(idx) mean(label[idx], na.rm = TRUE), numeric(1))
    Group_col <- avg_label
    color_title <- "Avg(label)"

  } else {
    lab_chr <- as.character(label)
    majority <- character(num_vertices)
    for (i in seq_len(num_vertices)) {
      pts <- piv[[i]]
      if (length(pts) > 0) {
        ux <- unique(lab_chr[pts])
        majority[i] <- ux[which.max(tabulate(match(lab_chr[pts], ux)))]
      } else {
        majority[i] <- "NA"
      }
    }
    Group_col <- factor(majority)
    color_title <- "Majority label"
  }

  if (!is.null(legend_name)) {
    color_title <- legend_name
  }

  # 3D force-directed layout of the Mapper graph. Edges with more shared
  # points get a higher weight, pulling their endpoints closer together --
  # the 3D analogue of forceNetwork's linkDistance/charge behaviour.
  set.seed(seed)
  if (num_vertices > 0) {
    fr_weights <- if (nrow(edge_list) > 0) edge_weights + 1 else NULL
    coords <- igraph::layout_with_fr(g, dim = 3, weights = fr_weights)
  } else {
    coords <- matrix(numeric(0), nrow = 0, ncol = 3)
  }

  node_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3],
    Nodename = seq_len(num_vertices),
    Nodesize = vertex.size * 5,
    Group = Group_col
  )
  node_df$marker_size <- pmax(6, sqrt(node_df$Nodesize) * 2)

  is_continuous <- is.numeric(node_df$Group)

  node_df$hover <- sprintf(
    "Vertex: %d<br>Points: %d<br>%s: %s",
    node_df$Nodename, vertex.size, color_title,
    if (is_continuous) round(node_df$Group, 3) else as.character(node_df$Group)
  )

  # NA-separated line segments let every edge be drawn in a single trace.
  if (nrow(edge_list) > 0) {
    edge_x <- as.vector(rbind(coords[edge_list[, 1], 1], coords[edge_list[, 2], 1], NA))
    edge_y <- as.vector(rbind(coords[edge_list[, 1], 2], coords[edge_list[, 2], 2], NA))
    edge_z <- as.vector(rbind(coords[edge_list[, 1], 3], coords[edge_list[, 2], 3], NA))
  } else {
    edge_x <- edge_y <- edge_z <- numeric(0)
  }

  axis_style <- list(
    title = "", showticklabels = FALSE, showgrid = FALSE,
    zeroline = FALSE, showbackground = FALSE
  )

  p <- plotly::plot_ly()

  if (length(edge_x) > 0) {
    p <- plotly::add_trace(
      p,
      x = edge_x, y = edge_y, z = edge_z,
      type = "scatter3d", mode = "lines",
      line = list(color = "rgba(120,120,120,0.6)", width = 2),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }

  if (is_continuous) {
    p <- plotly::add_trace(
      p,
      data = node_df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(
        size = ~marker_size,
        color = ~Group,
        colorscale = "Viridis",
        colorbar = list(title = color_title),
        line = list(color = "rgba(0,0,0,0.3)", width = 0.5)
      ),
      text = ~hover,
      hoverinfo = "text",
      showlegend = FALSE
    )
  } else {
    for (grp in levels(factor(node_df$Group))) {
      sub_df <- node_df[as.character(node_df$Group) == grp, , drop = FALSE]
      p <- plotly::add_trace(
        p,
        data = sub_df,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = grp,
        marker = list(
          size = ~marker_size,
          line = list(color = "rgba(0,0,0,0.3)", width = 0.5)
        ),
        text = ~hover,
        hoverinfo = "text",
        showlegend = legend
      )
    }
  }

  p <- plotly::layout(
    p,
    scene = list(xaxis = axis_style, yaxis = axis_style, zaxis = axis_style),
    showlegend = legend,
    legend = list(title = list(text = color_title))
  )

  return(p)
}

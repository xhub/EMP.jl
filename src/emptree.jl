# TODO: this should go to ReSHOP!

"""
Add a node of type MP


# Arguments
- emp: the EMP object
- mdl: an opaque pointer to the ReSHOP.reshop_model
- node: the node to process
- reshop_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMPmaster, mdl::Ptr{ReSHOP.reshop_model}, node::MathPrgm, reshop_mp_parent::Ptr{ReSHOP.mathprgm}=Ptr{ReSHOP.mathprgm}(C_NULL))
    reshop_mp = ReSHOP.reshop_declare_mathprgm(node, mdl, emp.NLoffset)
    node.solverobj = reshop_mp

#    ReSHOP.reshop_add_mp(reshop_emp, reshop_node, reshop_mp)

    if reshop_mp_parent != C_NULL
        ReSHOP.emp_add_mp_mp(reshop_mp_parent, reshop_mp)
    else
        ReSHOP.emp_set_root(mdl, reshop_mp)
    end

    for mp_mp in node.mps
        empnode(emp, mdl, mp_mp, reshop_mp)
    end

    for mp_equil in node.equils
        empnode(emp, mdl, mp_equil, reshop_mp)
    end

end

"""
Add a node of type equilibrium


# Arguments
- emp: the EMP object
- mdl: an opaque pointer to the ReSHOP.reshop_model
- node: the node to process
- reshop_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMPmaster, mdl::Ptr{ReSHOP.reshop_model}, node::Vector{MathPrgm}, reshop_mp_parent::Ptr{ReSHOP.mathprgm}=Ptr{ReSHOP.mathprgm}(C_NULL))
    reshop_mpe = ReSHOP.emp_create_equil(length(node))

    reshop_mps = []
    for mp in node
        reshop_mp = ReSHOP.reshop_declare_mathprgm(mp, mdl, emp.NLoffset)
        mp.solverobj = reshop_mp
        ReSHOP.emp_equil_add(reshop_mpe, reshop_mp)
        push!(reshop_mps, reshop_mp)
    end

    if reshop_mp_parent != C_NULL
        ReSHOP.emp_add_mp_equil(mdl, reshop_node, reshop_mpe)
    else
        ReSHOP.emp_set_root(mdl, reshop_mpe)
    end

    for (mp, reshop_mp) in zip(node, reshop_mps)
        for mp_mp in mp.mps
            empnode(emp, mdl, mp_mp, reshop_mp)
        end

        for mp_equil in mp.equils
            empnode(emp, mdl, mp_equil, reshop_mp)
        end
    end
end

"""
    emptree(emp, mdl)

Pass the EMP graph structure to the EMP solver
"""
function emptree(emp::EMPmaster, mdl::Ptr{ReSHOP.reshop_model})
    if !isnothing(emp.root)
        ReSHOP.emp_mp_ensure(mdl, length(emp.mps))

        empnode(emp, mdl, emp.root)
    end
end

import React from 'react';
import {NavLink, useNavigate} from 'react-router-dom';

export default function NavBar() {

    const navigate = useNavigate();
    const navigateToSignIn = () =>{
        //TODO change link to sign in page
        navigate('/');
    }

    return (
        <nav className='container nav-bar'>
            <h1>Bioinformatics Tools</h1>
            <ul>
                <li><NavLink to='/'>Home</NavLink></li>
                <li><NavLink to='/tools'>Tools</NavLink></li>
                <li><NavLink to='/about'>About</NavLink></li>
                <li><button className='sign-in-button' onClick={navigateToSignIn}>Sign In</button></li>
            </ul>
        </nav>
    )
}
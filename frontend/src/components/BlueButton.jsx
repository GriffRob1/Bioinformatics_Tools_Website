import React from 'react'
import {useState} from "react";
import {useNavigate} from "react-router-dom";

export default function BlueButton({content, URL, buttonClass}) {

    const [isMouseOver, setIsMouseover] = useState(false);

    const setBorderWhite = () => {
        setIsMouseover(true);
    }

    const setBorderBlack = () => {
        setIsMouseover(false);
    }

    const buttonStyle = {
        borderColor: isMouseOver ? 'white' : 'black'
    }

    const navigate = useNavigate();
    const navigateTo = () => {
        navigate(URL);
    }

    return (
        <button className={`blue-button ${buttonClass}`}
                onClick={navigateTo}
                onMouseOver={setBorderWhite}
                onMouseOut={setBorderBlack}
                style={buttonStyle}
        >
            {content}
        </button>
    )

}
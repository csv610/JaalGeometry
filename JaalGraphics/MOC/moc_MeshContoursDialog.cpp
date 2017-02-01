/****************************************************************************
** Meta object code from reading C++ file 'MeshContoursDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshContoursDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshContoursDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshContoursDialog_t {
    QByteArrayData data[11];
    char stringdata0[120];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshContoursDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshContoursDialog_t qt_meta_stringdata_JMeshContoursDialog = {
    {
        QT_MOC_LITERAL(0, 0, 19), // "JMeshContoursDialog"
        QT_MOC_LITERAL(1, 20, 14), // "displayContour"
        QT_MOC_LITERAL(2, 35, 0), // ""
        QT_MOC_LITERAL(3, 36, 11), // "getOriginal"
        QT_MOC_LITERAL(4, 48, 9), // "getSmooth"
        QT_MOC_LITERAL(5, 58, 10), // "getReparam"
        QT_MOC_LITERAL(6, 69, 11), // "getSimplify"
        QT_MOC_LITERAL(7, 81, 11), // "closeDialog"
        QT_MOC_LITERAL(8, 93, 13), // "keyPressEvent"
        QT_MOC_LITERAL(9, 107, 10), // "QKeyEvent*"
        QT_MOC_LITERAL(10, 118, 1) // "e"

    },
    "JMeshContoursDialog\0displayContour\0\0"
    "getOriginal\0getSmooth\0getReparam\0"
    "getSimplify\0closeDialog\0keyPressEvent\0"
    "QKeyEvent*\0e"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshContoursDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    7,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   49,    2, 0x08 /* Private */,
    3,    0,   50,    2, 0x08 /* Private */,
    4,    0,   51,    2, 0x08 /* Private */,
    5,    0,   52,    2, 0x08 /* Private */,
    6,    0,   53,    2, 0x08 /* Private */,
    7,    0,   54,    2, 0x08 /* Private */,
    8,    1,   55,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 9,   10,

    0        // eod
};

void JMeshContoursDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshContoursDialog *_t = static_cast<JMeshContoursDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->displayContour();
            break;
        case 1:
            _t->getOriginal();
            break;
        case 2:
            _t->getSmooth();
            break;
        case 3:
            _t->getReparam();
            break;
        case 4:
            _t->getSimplify();
            break;
        case 5:
            _t->closeDialog();
            break;
        case 6:
            _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1])));
            break;
        default:
            ;
        }
    }
}

const QMetaObject JMeshContoursDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshContoursDialog.data,
        qt_meta_data_JMeshContoursDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshContoursDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshContoursDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshContoursDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshContoursDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshContoursDialog"))
        return static_cast< Ui::MeshContoursDialog*>(const_cast< JMeshContoursDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshContoursDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
